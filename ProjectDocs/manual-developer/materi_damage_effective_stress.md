# materi_damage_effective_stress (damage stress fix, 2026-09-29)

Legacy defect of the 2014 sfnet source, found by the keyword convergence
cycle (the `materi_damage_apply` measurements of 2026-09-27 showed the
damaged stress differing from the Professional while `dam` was identical).

## The defect

`set_stress()` (stress.cc) evaluated the damage and scaled the stress in the
membrane loop:

```c
damage( gr, new_epe, new_sig, old_damage, new_damage );
array_multiply( new_sig, new_sig, 1.-new_damage, MDIM*MDIM );
```

At that point `new_sig` = `rotated_old_sig` + C·Δε, and `rotated_old_sig` is
the stress dof of the previous step, which is already the DAMAGED stress.
So each step multiplied the carried stress by (1−d) again:

    sigma_n = (1 - d_n) (sigma_{n-1} + C deps)

instead of the manual Professional 2.2.5 law `sigma_damaged = (1-d)
sigma_undamaged`. The damage compounded over the steps: the `damage1` bar
(E = 23.4e9, eps = 6e-4 in 60 steps) gave 1436.83 instead of
(1−d)·E·eps = 85598.32178. Holding the strain constant (Δε = 0, d constant)
made the stress decay to ~0. The recursion above reproduces 1436.828 exactly.
The same block is in `Sources-2014.zip`, so every GNU version had it.

## The fix: effective stress as integration-point state

The EFFECTIVE (undamaged) stress is kept as its own per-integration-point
state and evolved with the normal constitutive update. The damaged stress
is derived from it, `sigma = (1-d) sigma_eff`, and it is never divided by
(1−d). `materi_stress` keeps its meaning (the damaged stress, as in the
Professional).

### Storage

- Record `ELEMENT_INTPNT_MATERI_STRESS_EFFECTIVE`
  (`element_intpnt_materi_stress_effective`), registered in `database.cc`
  next to `element_intpnt_materi_undrained_pressure`, whose pattern it
  follows: class ELEMENT, DOUBLE_PRECISION, `version_all = 1` (VERSION_NORMAL
  = converged previous step, VERSION_NEW = current iterate; the step-end
  `db_version_copy( VERSION_NEW, VERSION_NORMAL )` in `top.cc` promotes it),
  `fixed_length = 0`, `data_length = MPOINT*MSTRESS_EFF_SLOTS`.
- **Internal**: `external = 0`, so it is neither read from the input file
  nor written by `print_database`. The Professional manual has no record for
  the undamaged stress (checked: no `element_intpnt_*` or other record for
  it), so no user keyword was invented.
- Layout: `MSTRESS_EFF_SLOTS` (7, `tochnog.h`/`tochnog-mod.h`) values per
  integration point: a flag (1 = stored) and the 6 stresses in
  `stress_indx()` order. The stress is stored in the same frame as the
  stress dofs (after the `memory` rotation).
- Enum appended at the end (`tochnog.h` and `tochnog-mod.h` in sync).
- Pre-allocation (zeros, flag 0) in `top.cc`: in `top()` for the initial
  elements and in `step_start()` every step for elements created later
  (macro meshes, refinement), only when `materi_damage` and `materi_stress`
  are initialized. `materi()` runs inside the parallel element loop and must
  not allocate; it does a read-modify-write of the VERSION_NEW record and
  PUTs it back with the existing length.

Why not the alternatives: the stress dof is a nodal (smoothed) quantity and
cannot carry a second stress without a user-visible dof; `hisv` belongs to
the user models; `old_sig/(1-d_old)` would divide by (1−d) (excluded by
design, singular for d → 1).

### Data flow

- `materi()` (materi.cc), after reading the old dofs: `old_sig_eff` = the
  stored effective stress of this integration point (VERSION_NORMAL, flag 1)
  or, without a stored value (first step, new element), the stress dof
  itself (d = 0 at the start). `new_sig` (the stress seen by stress-dependent
  material data in `set_stress`) is initialized with it. After the `memory`
  back-rotation `rotated_old_sig` is replaced by the same rotation of
  `old_sig_eff`. `old_sig` itself stays the damaged dof value: it feeds the
  stress-dof residual `(new_sig - old_sig)/dtime`, `materi_dynamic` and
  `materi_work`.
- `set_stress()` (stress.cc): new output argument `new_sig_eff[]`. In the
  membrane loop the elastic/viscous/... update runs on the effective stress;
  the damage block copies it to `new_sig_eff`, calls `damage()` with it and
  sets `new_sig = (1-d) new_sig_eff`; the membrane matrix `memmat` keeps its
  (1−d) scaling (it pairs with the damaged `new_sig` in `membrane_apply`).
  Without a damage evaluation (`damage_done` 0: no damage model, gated step,
  incremental plasticity branch) `new_sig_eff = new_sig`.
  `group_materi_factor` scales both.
- `materi()`: `new_sig_eff` gets the same forward rotation as `new_sig`
  and is written into the VERSION_NEW record (flag 1).

### Where plasticity acts

The plastic yield test stays on the DAMAGED stress (unchanged position of
the damage block, before the yield test). The standard effective-stress
formulation (yield on sigma_eff) was checked against the Professional
02-08-2026 and rejected: a 1D bar with Mazars damage and
`group_materi_plasti_vonmises 0 1.e7` (above the damaged peak stress,
below the effective stress reached) gives NO plastic strain in the
Professional (`sigxx` 85598.32, `eppxx` 0), which only the damaged-stress
yield test reproduces (GNU after the fix: identical). With
`vonmises 3.e6` (plasticity before damage) the Professional and the GNU
differ in the plastic strain path (`eppxx` 1.368e-4 vs 1.879e-4, `dam`
0.96946 vs 0.92537; the old GNU had the same strains as the new one). Each
binary satisfies `sigma = (1-d) C eps_e` with its own strains; the
remaining difference is not in the stress law (an iteration-order cause
was suspected first). Resolved 2026-10-12: it is the Professional equivalent
strain defect of BACKLOG C7 (class A; guard test
`test-2026/damage_mazars_vonmises1`).

### Damage switched off in a step

`materi_damage_apply` / `control_materi_damage_apply -no` (manual 6.799:
"any damage data ... will be ignored"): `set_stress` sets
`new_damage = old_damage` and skips the (1−d) scaling, and `materi()` does
not substitute the effective stress, so the step runs on the damaged stress
like a model without damage data; that stress is stored as the effective
one. Measured on the Professional 02-08-2026 with
`materi_damage_apply2.dat` (damage in timestep 0, gated in timestep 1):
`sigxx` 1.0872657685e7 in both (frozen-d semantics would give 7.705e6).
Re-enabling damage after a gated phase with d > 0 starts from that stored
stress (the Professional presumably divides by (1−d); not measured, not
reproduced).

## Verification

Pro = Professional 02-08-2026; theory = manual 2.2.5 with the Mazars d.

| model | theory | Pro | GNU old | GNU new |
|---|---|---|---|---|
| `damage_stress1` 1D, eps 6e-4, `-total_linear` | 85598.32178 | 85598.32178 | 1436.828 | 85598.32178 |
| `damage_stress_hold1` load + hold | 85598.32178 | 85598.32178 | -3.0e-13 | 85598.32178 |
| 1D `-updated`, no `materi_displacement`, load | (1−d)E·eps | 85598.32 | 1439.82 | 85800.49 (*) |
| same, load + hold | (1−d)E·eps | 85598.32 | -3.2e-14 | 85800.49 (*) |
| `damage_stress_3d1` hex8 load + hold, `sigyy` | — | -275960.2381 | -0.766 | -275960.2381 |
| same, load only, `sigyy` | — | -275960.2381 | -35511.71 | -275960.2381 |
| hex8 `diag(-1e-3,-1e-3,8e-4)`, `sigxx` | -673833.77 | -673833.77 | -88435.64 | -673833.77 |
| `materi_damage_apply2` (gated step 2) | — | 1.0872657685e7 | 284630.4 | 1.0872657685e7 |
| 1D + `vonmises 1.e7` | 85598.32 | 85598.32 | 1436.83 | 85598.32 |

(*) The `-updated` bar without `materi_displacement` reaches a total strain
5.99823e-4 in the GNU (6.0e-4 in the Professional), so d = 0.9938871 vs
0.9939033; the GNU stress equals (1−d)·E·eps for its own strain exactly
(85800.488). This strain-measure difference is pre-existing (identical in
the old binary) and independent of this fix.

Self-checking inputs (`validation-suite/test-2026/`): `damage_stress1.dat`,
`damage_stress_hold1.dat`, `damage_stress_3d1.dat` (pass on the GNU and the
Professional; fail with rc=1 on the pre-fix binary and when a target is
changed).

Regression (clean build): validation-suite 233 PASS / 1 RUNFAIL and sfnet
187 PASS / 10 RUNFAIL, classification identical to the pre-fix binary;
test-2026 15/15; `build_safe.sh` 16/16 runs, same 70 file-check failures.
The only real `.dbs` change in both corpora is sfnet `damage1` (`sigxx`
1436.83 -> 85598.32; its target checks `dam` only, unchanged). The other
`.dbs` differences are uninitialized `element_dof` padding (denormal
garbage) and time stamps, identical in nature between two runs.
Existing targets: `materi_damage_apply1/2.dat` check `dam` (and `sigxx` with
damage off), `materi_failure_apply*.dat` have no damage: no target encoded
the buggy stress.

## Mazars gamma (checked, not changed)

`damage_mazars()` computes the triaxiality factor from `new_sig[0..2]`,
i.e. xx, xy, xz of the 3x3 matrix, not the principal stresses. It now
receives the effective stress (gamma is scale-invariant, so proportional
paths give the same d as before). Probed on single hex8 elements:

| strain state | GNU d | Pro d | GNU gamma form reproduces GNU |
|---|---|---|---|
| diag(-1e-3,-1e-3,8e-4) | 0.4224282004 | 0.4224282004 | yes (principal gamma would give 0.2385) |
| yy -5e-4, zz 3e-4, yz 5e-4 | 0.4156136134 | 0.4156136134 | yes |
| pure shear xy 1e-3 | 0.6524270351 | 0.6524270351 | yes |
| xx -5e-4, xy -5e-4, zz 3e-4 | 0.1080924614 | 0.1875189548 | yes |
| xx -5e-4, xy +5e-4, zz 3e-4 | 0.2635138083 | 0.3935182025 | yes |

The GNU form is not objective (the last two states are reflections of each
other and give different d), but the Professional is not objective either
and matches neither the GNU form, nor the principal-stress gamma, nor
gamma = 1 (the Professional manual 2.2.5 has no gamma at all). A fitted
equivalent gamma is 0.8634 and 1.2187 for the two states. With no evidence
for a specific correct form, gamma was left unchanged; it only matters
when a compressive xx is combined with a same-sign xy or xz. Open item in
`BACKLOG-OPEN-ISSUES.md`.

Resolved 2026-10-12 (F4 T6f, BACKLOG C7): the Professional difference in
the last two states is not gamma but the equivalent strain. Both states
have two positive principal strains, and the Professional accumulates
them as s = 2 s + eps^2 in ascending order (weights 4, 2, 1; the 2014
`epseq += epseq + ...` defect fixed in the GNU on 2026-08-10, see
`group_materi_damage_mazars.md`), which reproduces 14 probes to 7e-14
with the GNU gamma. The plastic-path difference (`vonmises 3.e6`) is the
same defect. Deliberate divergence, class A (owner 2026-10-12).

## Known limitations / pending

- Behaviour change beyond the defect: for damage models the step update
  now starts from the integration-point effective state instead of the
  nodal stress dofs interpolated to the integration point, so in
  non-homogeneous meshes the stress history is no longer smoothed through
  the nodes (homogeneous single-element models: identical).
- Direct stress cut-offs (`materi_direct_full_mc`, `materi_direct_cutoff`,
  `materi_compression_cutoff`) act on the damaged stress after the membrane
  loop; the effective stress does not see the cut (it would need a
  division). Damage combined with these cut-offs is not verified.
- The effective state is not in the printed database, so a restart from a
  printed database starts from the damaged stress dofs (like the first
  step).
- `ddsdde` (the momentum tangent) is still the undamaged `Cmem` (legacy,
  unchanged): convergence rate only, not the converged stress.
- Remeshing/renumbering move element records; an element whose record is
  missing restarts from its stress dofs.

## Differences with the Professional version

- **Mazars equivalent strain** (deliberate divergence, class A, project
  decision 2026-10-12; BACKLOG C7). GNU (`damage_mazars()`, `damage.cc`): manual 2.2.5, eps_eq =
  sqrt(sum <eps_i>+^2) of the positive principal strains (hex8
  diag(2e-4,-3e-4,1e-4): dam 0.28610207515). Professional 02-08-2026:
  accumulates them in ascending order as s = 2 s + eps^2 (weights 4, 2, 1;
  the 2014 accumulation defect the GNU fixed on 2026-08-10): dam 0.36409.
  It differs whenever two or three principal strains are positive, also
  through a von Mises + damage path (bar with `vonmises 3e6`: GNU dam
  0.925372, eppxx 1.87863e-4; Professional 0.96946, 1.368e-4). Guard
  tests `test-2026/damage_mazars_eq_two_positive1`,
  `damage_mazars_vonmises1`. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- The triaxiality factor gamma is built from the xx, xy, xz stress
  components in both codes (`damage.cc:84-95`) (not objective: a compressive xx with xy or xz
  shear gives an orientation-dependent `dam`). Same in the Professional,
  so not a divergence; noted in BACKLOG C7, no fix planned.
- **`-updated` without `materi_displacement`: total strain (damage bar)**
  (known difference, not decided; AUDIT section 5 (known items)). GNU:
  5.99823e-4. Professional 02-08-2026: 6.0e-4 (dam and stress differ ~0.2
  %). Status: Pre-existing, not decided. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
