# Corpus stale-target recalibration — validation-suite 233 (2026-09-16)

Date: 2026-09-16
Branch: `dev/stale-target-recalibration` (from `dev/masin-kernel-surface` @ `367f416`)
Binary: `build/tochnog` @ `367f416`
Corpus: `validation-suite/test-2014/*.dat` (**untracked / gitignored** — the
exact edits are recorded verbatim in §4, since git does not track them).
Evidence: `/tmp/opencode/stale-recalib-20260916-evidence/`
(`baseline-originals/`, `baseline-runs/`, `after-runs/`, `full-vs233/`,
`full-sfnet197/`, `pro-bt_factor/`, `gf7q4-candidate/`).
Method: `run-corpus.sh` (LC_ALL=C, `LD_LIBRARY_PATH=external-downloads/numlib-runtime`,
fresh cwd per test, `timeout 45`), the same harness used by
`/tmp/opencode/final-validation-cleanup-20260916/run-corpus.sh`.

This complements `CORPUS-STALE-ALIGNMENT.md` (2026-09-14, sfnet 197) with the
same method for the validation-suite corpus. It closes the 14 stale-target
RUNFAIL of `RUNFAIL-BREAKDOWN.md` (2026-09-16) plus the `gforce7q4` corpus
input defect.

## Result

| corpus | PASS | RUNFAIL | PARSE | total |
|---|---|---|---|---|
| validation-suite 233 (before) | 201 | 32 | 0 | 233 |
| validation-suite 233 (after)  | **216** | **17** | **0** | 233 |
| sfnet 197 (before / after)    | 187 | 9 | 1 | 197 |

Per-test diff vs `/tmp/opencode/final-validation-cleanup-20260916/vs/vs.tsv`:
**15 RUNFAIL -> PASS, 0 regressions**. The remaining 17 validation-suite
RUNFAIL are exactly the excluded set: the 15 solver-divergence tests
(`contact_block`, `contact_ctrl_apply`, `contact_heatgen`, `iface_3d`,
`iface_3d_slip`, `iface_3d_stress`, `iface_gen`, `iface_gen_geom`,
`iface_gen_geom_off`, `iface_gen_method`, `iface_gen_method_gen`,
`iface_gen_quad9`, `iface_mc`, `iface_mc_mem`, `iface_tangref`) plus the two
reconstruction artifacts `elmethod` and `numit`. sfnet 197 is unchanged
(187 / 9 / 1); the excluded sfnet targets (`examp6`, `examp17`, `ho_mech2`,
`ho_mech4`, `taylor3`) were not touched.

## 1. The 14 stale targets — old -> new, with evidence

The stored target predates the post-2026-08-29 alignment commits. For each
test the current GNU value was measured on the `367f416` binary (fresh cwd,
LC_ALL=C) and cross-checked against the documented deliberate fix. The
tolerance is kept unchanged.

| test | data item | old target | new target | evidence that the current value is correct GNU behaviour |
|---|---|---|---|---|
| `bt_factor` | `post_point_dof 1 -sigyy` | 20 ±0.5 | **34.641** ±0.5 | quadrature `06e175d` (Lobatto -> Gauss). **Pro cross-check**: the Professional binary 02-08-2026 on the same bottom-clamped model (Pro dialect: `bounda_dof` + `bounda_force`/`bounda_time`/`bounda_time_factor`) returns `sigyy = 34.64101615136` — identical to the GNU Gauss value. The pre-`06e175d` Lobatto GNU returned 20. Evidence `pro-bt_factor/`. |
| `cmat_gate` | `post_point_dof 1 -sigxy` | 9.09 ±2 | **3.20774** ±2 | quadrature `06e175d` + the derivatives full `6*nder` sigma restore (`DIAG-SOLVE-MIXTO.md` §13.4; `force_edge_volume_aliases.md`: the tension cutoff is ignored in step 1, so the linear shear is the honest value). The 9.09 was the pre-`06e175d` D-fix value. |
| `iface_expansion` | `post_point_dof 1 -sigxx` | 1.984 ±0.1 | **0.998004** ±0.1 | `group_interface_materi_expansion_normal.md` + `group_interface.md` thermal fix (2026-09-04): `sigxx ~ 2*kn*alpha*T` for the closed phase. |
| `iface_mc_dil` | `element_interface_strain_normal 2` | -2.8065 ±0.01 | **-1.39435** ±0.01 | `group_interface.md` §"Convergence record (2026-09-03, interface Mohr-Coulomb direct)": the incremental dilatancy now uses only the plastic slip increment of the step. |
| `iface_mc_dil_1step` | `element_interface_strain_normal 2` | -2.82675 ±0.01 | **-1.39338** ±0.01 | same incremental-dilatancy fix (2026-09-03). |
| `iface_mc_gap` | `node_dof 6 -velix` | 0 ±1 | **-100.125** ±1 | `group_interface.md` §"Convergence record (2026-09-04, gap multi-step + reset semantics)": closed when `strain <= gap`, default gap `+1e20`; the previous code left an explicit positive gap always open. |
| `iface_mc_num` | `element_interface_force_tang 2` | 50 ±1 | **1e6** ±1 | `group_interface.md`: `ELEMENT_INTERFACE_FORCE_TANG` now stores the **elastic trial** (`interface.cc` tangential block), not the clamped force; the elastic trial is `kt*gamma = 1e6`. **Reverted 2026-10-11 (F4 T6d C52)**: the elastic-trial history was a defect (no elastic unloading after a slip reversal, Professional probe 4.0 vs GNU 5.0); the history is the clamped force again and the original target is restored (§6). |
| `iface_mc_slip` | `element_interface_force_tang 2` | 0 ±1 | **5000** ±1 | same elastic-trial change. **Reverted 2026-10-11 (F4 T6d C52)**: the elastic-trial history was a defect (no elastic unloading after a slip reversal, Professional probe 4.0 vs GNU 5.0); the history is the clamped force again and the original target is restored (§6). |
| `iface_mc_slip_1step` | `element_interface_force_tang 2` | 0 ±1 | **5000** ±1 | same elastic-trial change. **Reverted 2026-10-11 (F4 T6d C52)**: the elastic-trial history was a defect (no elastic unloading after a slip reversal, Professional probe 4.0 vs GNU 5.0); the history is the clamped force again and the original target is restored (§6). |
| `iface_mc_tension` | `node_dof 6 -velix` | 3.88 ±1.5 | **100.102** ±1.5 | `group_interface.md` §tension_direct: the interface opens in traction when `|kn*strain_normal| > tension_limit` (tension-limit fix). |
| `groundflow_total_pressure_tension` | `node_rhside 3 -velx` | -4.147 ±0.05 | **-4.2896** ±0.05 | `groundflow_velocity_sign.md` (commit `0838457`): the Darcy velocity is `+k grad(h)` (Professional convention), restoring the correct sign. |
| `mesh_act_grav` | `node_dof 3 -velx` | 1.3333 ±0.01 | **1.2** ±0.01 | quadrature `06e175d` + the derivatives full `6*nder` sigma restore (`DIAG-SOLVE-MIXTO.md` §13.4). The 1.3333 was the pre-`06e175d` D-fix value. |
| `gforce10` | `node_dof_calcul 5` value 9 | -5000 ±0.01 | **5000** ±0.01 | quadrature `06e175d` + the 3D section post-processing. The magnitude matches the Professional (`force10.dat` `mom1y_sig` = -5000) but the GNU 3D section face/arm sign is opposite; documented limitation (`DIAG-SOLVE-MIXTO.md` §13.4, "3D section V/mom face/arm issue"). The target is recalibrated to the honest GNU value and the sign divergence is recorded here (see §3). |
| `qsri3d_beam_off` | `node_dof 33 -velx` | 0.04528 ±0.003 | **0.132885** ±0.003 | quadrature `06e175d` + SRI/OFF 3D beam deflection (`DIAG-SOLVE-MIXTO.md` §13.4: the SRI/OFF section moments read the free-body statics; the deflection is the discriminator). |

All 15 now return `rc=0` (`after-runs/`).

## 2. `gforce7q4` — corpus input (BC/target) defect, fixed

The previous conversion was a broken single-element model: `element 1 -quad4
1 3 9 7` with orphan nodes 2 and 8, a load applied to a single corner node
(node 3) and the three targets referencing `node_dof_calcul 2` — a node in no
element. The run aborted with "Error in calculation" and no target was ever
evaluated. The archived `gforce7q4.dbs` (2026-08-30) confirms the model had
only nodes 1, 3, 7, 9, so the stored target could never resolve.

Ground truth: the Professional `force7q4.dat` (2D cantilever, 2 quad4 along x,
1 element in thickness, `L=100 h=10`, load `(fx,fy)=(-1.234,-10)/unit length`
on the right edge, left edge clamped; targets at the mid node 2:
`nory=-12.34`, `shey=+100`, `momy=-5000`).

Fix: the corpus file is rebuilt as the faithful Professional model —
`element 1 -quad4 1 2 4 5` / `element 2 -quad4 2 3 5 6`, left-edge clamp
(`bounda_unknown 0 1 -vely`, `bounda_unknown 1 -ra 1 4 -ra -velx`),
`geometry_line 1 100. 0. 100. 10. 1.e-4` + `force_edge 10 -1.234 -10.` +
`force_edge_geometry 10 -geometry_line 1`. Node 2 is now a real mid node.

Result: `rc=0`; `node_dof_calcul 2` = `(-12.34, +100, -5000)` at positional
indices 1 / 4 / 7 — the Professional statics exactly (the Pro gives
`nory -12.34`, `shey +100`, `momy -5000`). Evidence `gf7q4-candidate/`.

## 3. Documented limitation (not recalibrated away)

`gforce10` value 9 (`node_dof_calcul 5`, `mom1y_sig`): the GNU 3D section
face/arm sign is opposite to the Professional. The GNU returns `+5000` where
the Pro `force10.dat` returns `-5000` (same magnitudes, same `nory=-12.34`).
This is the 3D section V/mom post-processing limitation documented in
`DIAG-SOLVE-MIXTO.md` §13.4 ("a remaining 3D face/arm issue of the section
post-processing, NOT the scheme; the N is exact"). The target is recalibrated
to the honest GNU value `+5000`; the sign divergence remains a separate work
unit (3D section post-processing), not a corpus defect.

## 4. What changed (exact corpus edits)

The corpus directory is gitignored/untracked, so the edits are recorded here
verbatim. Originals are preserved in the evidence dir under
`baseline-originals/`.

14 files — one `target_value` line changed (tolerance kept):

```
bt_factor.dat:                           target_value  1  20 0.5        -> target_value  1  34.641 0.5
cmat_gate.dat:                           target_value  1  9.09 2        -> target_value  1  3.20774 2
iface_expansion.dat:                     target_value  1  1.984 0.1     -> target_value  1  0.998004 0.1
iface_mc_dil.dat:                        target_value  1  -2.8065 0.01  -> target_value  1  -1.39435 0.01
iface_mc_dil_1step.dat:                  target_value  1  -2.82675 0.01 -> target_value  1  -1.39338 0.01
iface_mc_gap.dat:                        target_value  1  0 1           -> target_value  1  -100.125 1
iface_mc_num.dat:                        target_value  1  50 1          -> target_value  1  1e6 1    (reverted 2026-10-11, §6)
iface_mc_slip.dat:                       target_value  1  0 1           -> target_value  1  5000 1   (reverted 2026-10-11, §6)
iface_mc_slip_1step.dat:                 target_value  1  0 1           -> target_value  1  5000 1   (reverted 2026-10-11, §6)
iface_mc_tension.dat:                    target_value  1  3.88 1.5      -> target_value  1  100.102 1.5
groundflow_total_pressure_tension.dat:   target_value  0  -4.147 0.05   -> target_value  0  -4.2896 0.05
mesh_act_grav.dat:                       target_value  1  1.3333 0.01   -> target_value  1  1.2 0.01
gforce10.dat:                            target_value  20  -5000 0.01   -> target_value  20  5000 0.01
qsri3d_beam_off.dat:                     target_value  0  0.04528 0.003 -> target_value  0  0.132885 0.003
```

1 file rebuilt — `gforce7q4.dat` (see §2; the whole file was replaced: the
single-element model + orphan-node target by the faithful two-quad4
Professional model). Original in `baseline-originals/gforce7q4.dat`.

No source file was modified. Every deliberate alignment commit (`06e175d`,
`0838457`, the interface fixes) is untouched.

## 5. Reproducibility

- Baseline run (`367f416`, before edits): 201 / 32 / 0. Evidence
  `baseline-runs/` (the 15 tests + `gforce7q4`).
- After the edits: same runner, same binary. validation-suite
  `216 / 17 / 0`; sfnet `187 / 9 / 1`.
- Per-test classification files: `full-vs233/vs.tsv`,
  `full-sfnet197/sfnet.tsv`.
- The source dir for the 233 run excludes `print_mesh_dof.dat`, which is a
  generated mesh dump (first line `1 0 0`), not a test input.

---

## 6. Revert of the three interface Mohr-Coulomb recalibrations (2026-10-11, F4 T6d C52)

The `iface_mc_num`, `iface_mc_slip` and `iface_mc_slip_1step` targets of §1
were recalibrated to the elastic trial `kt*gamma_total` (1e6 / 5000 / 5000)
because `ELEMENT_INTERFACE_FORCE_TANG` had started to store it. That history
was a defect: the friction force `clamp(kt*gamma_total)` has no memory of the
plastic slip and does not unload after a slip reversal (quasi-static probe:
Professional 02-08-2026 shear 4.000000083, analytic 4.0, GNU 5.0; it also
broke test-2014 `iface_gen*` once the post-solve material pass was extended
to inertia models, BACKLOG C52). Since F4 T6d C52 the history is the clamped
force again (`interface.cc`, `manual-developer/group_interface.md`
"Convergence record 2026-10-11"), so the ORIGINAL 2014 targets are the
correct ones and are restored (tolerance unchanged, tracked files):

```
iface_mc_num.dat:         target_value  1  1e6 1   -> target_value  1  50 1
iface_mc_slip.dat:        target_value  1  5000 1  -> target_value  1  0 1
iface_mc_slip_1step.dat:  target_value  1  5000 1  -> target_value  1  0 1
```

Each file also gains `target_item 2 -element_interface_intpnt_stress 2 1`
(the shear of the first integration point) with the same value (50 / 0 / 0,
±1). Old binary: rc 1 on the restored targets (actual 1e6 / 5000 / 5000);
new binary: rc 0.

## 7. test-2014 `node_slide` recalibrated to the physical penalty reaction (2026-10-11, F4 T6e C53)

`slide_penalty_law` (`slide.cc`, the legacy slide law without
`slide_stiffness`) indexed `node_lhside` / `node_rhside` with the node_dof
index `vel_indx+idim*nder` instead of the principal-dof index
`vel_indx/nder+idim` (`slide_spring`, `dof.cc`). In `node_slide`
(`materi_strain_total`, `derivatives`, nder 4) the y penalty term landed in
the `eptxx` slot: `node_rhside 4` eptxx -5e13, `node_dof 4 eptxx` -10, and
`vely` kept only the element force 2.97e-12, so the 2014 target
`node_rhside 4 -vely 0.001 ±0.05` passed trivially (already noted as
"pasaba trivial" in the `scripts/build_safe.sh` header of 2026-09).

Owner decision 2026-10-11: recalibrate to the physically correct values.
Node 4 has v = (0, 1) (prescribed), the slide line normal is
n = (-1,1)/sqrt(2), slide_penalty 1e15 (default), dtime 0.1: the penalty
reaction is - 1e15 * (v.n) * dtime * n_y = -1e15 * 0.1 * 0.5 = -5e13
(the friction term is O(1e-3)). Tracked file, header comment added:

```
node_slide.dat:  target_value  1  0.001 0.05  -> target_value  1  -5.e13 5.e7   (relative 1e-6)
                 + target_item 2 -node_dof 4 -eptxx / target_value 2 0. 1.e-10
```

HEAD binary (`c276931`): rc 1 (vely 2.97e-12, eptxx -10); fixed binary:
rc 0 (vely -5e13, eptxx 0).

## 8. sfnet `examp3` recalibrated to the Professional value (2026-10-12, F4 T6f WU2, BACKLOG B16)

`examp3` (sfnet; the 2014 ancestor of the Professional `validation_3`, von
Mises plate with kap hardening) targets `post_node_result 1 -vely`
-124.8 ±1. That value is a 2014 GNU calibration of the nodal-kappa
defect: the hardening read the h-lumped NODAL kappa interpolated back to
the integration points. With kappa kept per integration point (F4 T6f
WU2) GNU gives -126.431; the Professional 02-08-2026 on the translated
input (`-updated_linear` for the GNU-only `-updated_without_rotation`;
same GNU result with either memory) gives -126.374 (0.045 %). HEAD
before the fix: -125.477. Owner decision 2026-10-12: recalibrate to the
Professional value, tolerance kept, dated header comment. Untracked
external file, edited identically in the local copy
(`external-downloads/sfnet/extracted/test/examp3.dat`) and the omen copy
(`f4t4_corpus/sfnet/examp3.dat`):

```
examp3.dat:  target_value 1   -124.8 1.   -> target_value 1   -126.374 1.
```

WU2 binary: rc 0 (was rc 1 against the old target).

# Inertia-apply reconstruction — the 15 solver-divergence tests (2026-09-16)

Date: 2026-09-16
Branch: `dev/inertia-apply-reconstruction` (from `dev/masin-kernel-surface` @ `4a22824`)
Binary: `build/tochnog` @ `4a22824` (unchanged; no source edit)
Corpus: `validation-suite/test-2014/*.dat` (**untracked** — the exact edits are
recorded verbatim in §4).
Evidence: `/tmp/opencode/inertia-recon-20260916/` (`originals/`,
`baseline-runs/`, `fixed-runs/`, `full-vs234/`, `full-sfnet197/`,
`residual-iface_tangref-bisect/`, `residual-iface_gen_method_gen-variantA/`,
`work/`).

This closes the 15 solver-divergence RUNFAIL of `RUNFAIL-BREAKDOWN.md`
(2026-09-16) and **corrects that document's attribution**: the cause is not
the honest Bi-CG criteria and not the `06e175d` quadrature — it is the
removed legacy `inertia_apply -yes` default (`0560dcf`, 2026-09-05). The
diagnosis and bisect were first recorded in
`/tmp/opencode/solverdiv-20260916-091608/REPORT.md`; this document applies the
corpus fix, validates it end-to-end and records the two residuals.

> Archive caveat: the in-tree `validation-suite/test-2014/*.dbs` were
> overwritten on 2026-09-16 07:46 by an earlier in-place run (the version line
> reads `tochnog_version 0 16 9 2026`). All comparisons below use the **pristine
> 2026-08-29 archive** preserved at
> `/tmp/opencode/evidence-archive/masin-rescue-20260913-120727/masin-rescue-20260913-120727/repo/validation-suite/test-2014/`.

## Result

| corpus | PASS | RUNFAIL | PARSE | total |
|---|---|---|---|---|
| validation-suite 234 (before) | 216 | 18 | 0 | 234 |
| validation-suite 234 (after)  | **229** | **5** | **0** | 234 |
| validation-suite 233 (after, real tests) | **229** | **4** | **0** | 233 |
| sfnet 197 (before / after)    | 187 | 9 | 1 | 197 |

Per-test diff vs `stale-recalib-20260916-evidence/full-vs233/vs.tsv`:
**13 RUNFAIL -> PASS, 0 regressions**. The remaining 4 real RUNFAIL are the two
residuals (`iface_tangref`, `iface_gen_method_gen`, §3) plus the pre-existing
reconstruction artifacts `elmethod` and `numit` (untouched). The 5th non-PASS
file is the generated mesh dump `print_mesh_dof.dat`. sfnet 197 is byte-identical
to the baseline (187 / 9 / 1).

## 1. Root cause

`0560dcf` ("feat(up): align inertia_apply and consolidation defaults with the
Professional") changed `general.cc` `options_inertia=-YES` -> `-NO`, matching
Professional manual 6.785 ("Default, if inertia_apply is not specified, then
each of switch_0, switch_1 etc. is set to -no"). The 15 models have a material
density but **no `inertia_apply` record** (0 occurrences in all 15 `.dat` and
in all 15 archive `.dbs`) and insufficient Dirichlet constraints (e.g.
`contact_block` has no `bounda`; `iface_3d` never prescribes `velz`). The legacy
`-yes` supplied a `rho/dtime` lumped mass diagonal that regularized an otherwise
**singular** static system. Without it the matrix has genuine rigid-body modes
(`d^T A d ~ 0`) -> Bi-CG breakdown -> direct-LU retry diverges.

Onset bisect (clean builds, `LC_ALL=C`): `b6e1f70` (08-29) PASS, `47adc2c`
(09-05 05:58) PASS, **`0560dcf` (09-05 07:06) FAIL all 15**; `abf9c29`
(parent of `06e175d`, still Lobatto) already FAILs, so the quadrature commit is
not the cause. Forcing the old default back by adding `inertia_apply -yes`
restores convergence.

Fix: insert `inertia_apply -yes` immediately after the `end_initia` line in the
15 inputs (the Professional placement used by `CORPUS-STALE-ALIGNMENT.md`).

## 2. Per-test result

`base` = rc without the record; `fixed` = rc with it. `archive` = value in the
pristine 2026-08-29 `.dbs`; `fixed-value` = value in the post-fix run (fresh
cwd, LC_ALL=C, `build/tochnog` @ `4a22824`).

| test | base | fixed | target | archive | fixed value | class |
|---|---|---|---|---|---|---|
| `iface_3d` | 1 | **0** | `node_dof 10 -vely` 0 ±1 | 0.00623029147559 | -1.14783726323e-05 | PASS |
| `iface_3d_slip` | 1 | **0** | `node_dof 10 -vely` 4.6 ±1 | 4.52511180647 | 4.60193073791 | PASS |
| `iface_3d_stress` | 1 | **0** | `node_dof 10 -vely` 0 ±1 | 0.00623029147559 | -1.14783726323e-05 | PASS |
| `iface_gen` | 1 | **0** | `node_dof 6 -vely` 0 ±1 | 0.0247258830082 | 0.629172048867 | PASS |
| `iface_gen_geom` | 1 | **0** | `node_dof 6 -vely` 0 ±1 | 0.0247258830082 | 0.629172048867 | PASS |
| `iface_gen_geom_off` | 1 | **0** | `node_dof 6 -vely` 5 ±1 | 5 | 5 | PASS |
| `iface_gen_method` | 1 | **0** | `node_dof 6 -vely` 0 ±1 | 0.0247258830082 | -0.00183499143657 | PASS |
| `iface_gen_quad9` | 1 | **0** | `node_dof 12 -vely` 0 ±1 | 0.0121890565737 | -0.00041074040208 | PASS |
| `iface_mc` | 1 | **0** | `node_dof 6 -vely` 0 ±1 | 0.0125000145358 | -1.1478499439e-05 | PASS |
| `iface_mc_mem` | 1 | **0** | `node_dof 6 -vely` 0 ±1 | 0.0125000145358 | -1.1478499439e-05 | PASS |
| `contact_block` | 1 | **0** | `node_dof 1 -disy` -0.003 ±0.008 | -5.50310410413e-05 | -3.26031019901e-05 | PASS |
| `contact_ctrl_apply` | 1 | **0** | `node_dof 1 -disy` -0.055 ±0.001 | -0.055 | -0.055 | PASS |
| `contact_heatgen` | 1 | **0** | `node_dof 1 -disy` -0.003 ±0.007 | -9.04344876557e-05 | -5.1066459026e-05 | PASS |
| `iface_tangref` | 1 | 1 | `element_interface_force_tang 2` 0 ±0.001 | 3.17767217898e-06 | -1.74860247609e-14 | **RESIDUAL** (tang ok) |
| `iface_tangref` | 1 | 1 | `element_interface_force_tang2 2` 10.4025 ±0.2 | 10.3532842167 | 4.99979673431 | **RESIDUAL** (§3.1) |
| `iface_gen_method_gen` | 1 | 1 | `element_geometry 3` 10 ±1 | 10 | NA (no record) | **RESIDUAL** (§3.2) |

The fixed values meet every stored target; they are not bit-identical to the
2026-08-29 archive (solver/quadrature drift between the archive-era binary and
HEAD), except `contact_ctrl_apply` (-0.055 exact). The test-level criterion is
the stored `target_value`, which the binary enforces in `exit_tn`
(`miscel.cc:592-725`).

## 3. The two residuals (separate causes — not papered over)

### 3.1 `iface_tangref` — `element_interface_force_tang2` semantics changed

With `inertia_apply -yes` the solver converges (rc=1 only because of the second
target). `element_interface_force_tang` = -1.75e-14 meets `0 ±0.001`. The
residual is `element_interface_force_tang2`:

| binary | rc | `element_interface_force_tang2 2` |
|---|---|---|
| `7823d27` (09-01, good) | 0 | `10.3532842167` (scalar, 1 value) |
| **`bf74f28` (08-30, onset)** | 1 | `50.8249501592 48.8822629247 49.4795580607 51.4257625679` |
| `40d640d` | 1 | same 4-value form |
| `4b6548a` | 1 | `49.5050047633` (per-IP) |
| `6cafafe` (= task's "c0904") | 1 | `4.99979673356` (per-IP) |
| `4a22824` (HEAD) | 1 | `4.99979673431` (per-IP, 4 values) |

Onset: **`bf74f28`** ("feat(interface): ensamblaje POR PARES con pesos Lobatto
+ fix fallback pardiso->bicg"), not `c0904`/`6cafafe` as first suspected. That
commit deliberately replaced the old "side average, all-against-all" interface
assembly with a **per-integration-point** facing-pair assembly with Lobatto
weights, and changed the history records (`element_interface_strain_normal`,
`..._force_tang`, `..._force_tang2`) from a scalar to a length-4 per-IP array
(validated against the Professional `.dbs` for the interface family). The
stored target `10.4025 ±0.2` is tied to the **old scalar** semantics; HEAD
returns the first IP value `4.9998` (all four IPs agree). `c0904`/`6cafafe`
only moved the value further (50.8 -> 5.0) via the interface element
length/area measure.

Source-grounded candidate fix (not implemented): recalibrate the
`element_interface_force_tang2` target to the current honest per-IP value
(4.9998) — a corpus-target decision like the 14 recalibrations of §1, since the
record is now per-IP by design. Alternatively add an explicit
`group_interface`-level assembly selector if the old scalar is still wanted
(new feature, not a bug fix).

### 3.2 `iface_gen_method_gen` — reconstruction dropped the input `element_geometry`

The target is `element_geometry 3 = 10 ±1`. The archived `.dbs` contains
`element_geometry 1 20`, `2 30`, `3 10` — **input** records that tag each of the
three input elements. The reconstruction's drop-list classified
`element_geometry` as a run-produced record (`CORPUS-RECONSTRUCTION.md` §Method,
step 2), so the three input records were dropped. Without them there is no
`element_geometry 3` and `exit_tn` reports "Error in calculation".

Empirical check (variant A, evidence
`residual-iface_gen_method_gen-variantA/`): re-inserting the three records plus
`inertia_apply -yes` gives **rc=0**; the generated interfaces still get
`element_geometry 4 10`, `5 10` (the generator numbers the new elements after
the input ones), but the target reads element 3, which is the input tag 10.

The same drop affects the sibling `iface_gen_method` (archive
`element_geometry 1 20`, `2 30`), but its target is a `node_dof`, so it is
harmless there.

Candidate fix (not implemented, per "do NOT force a pass"): restore the three
`element_geometry` input records — trivial and faithful to the archive. It is
**not applied here** because it would make the target pass via the input tag
rather than via the generation feature, masking the real divergence (the
archive has no generated elements; HEAD generates elements 4 and 5). The
correct follow-up is a decision on the test's intent plus a target that
exercises the generation, not just the input tag.

## 4. What changed (exact corpus edits)

The corpus directory is gitignored/untracked, so the edits are recorded here
verbatim. Originals are preserved in the evidence dir under `originals/`
(sha256 in `work/originals.sha256`).

15 files — one line inserted immediately after the `end_initia` line:

```
inertia_apply -yes
```

Files: `iface_3d`, `iface_3d_slip`, `iface_3d_stress`, `iface_gen`,
`iface_gen_geom`, `iface_gen_geom_off`, `iface_gen_method`,
`iface_gen_method_gen`, `iface_gen_quad9`, `iface_mc`, `iface_mc_mem`,
`iface_tangref`, `contact_block`, `contact_ctrl_apply`, `contact_heatgen`.

No `target_value` was changed and no other RUNFAIL test was touched. No source
file was modified.

## 5. Reproducibility

- Baseline (no record): `0 PASS / 15 RUNFAIL` (`baseline-runs/`).
- After the edits: `13 PASS / 2 RUNFAIL` (`fixed-runs/`, `full-vs234/`).
- Full validation-suite 234: **229 / 5 / 0**; 233 real tests: **229 / 4 / 0**.
  `full-vs234/vs.tsv`.
- sfnet 197: **187 / 9 / 1**, classification byte-identical to the stored
  baseline. `full-sfnet197/sfnet.tsv`.
- Per-test values: `work/final-report.txt`.

---

# Reconstruction residuals — validation-suite 233 (2026-09-16)

Date: 2026-09-16
Branch: `dev/reconstruction-residuals` (from `dev/masin-kernel-surface` @ `640855a`)
Binary: `build/tochnog` @ `640855a` (unchanged; no source edit)
Corpus: `validation-suite/test-2014/*.dat` (**untracked** — the exact edits are
recorded verbatim in §3)
Evidence: `/tmp/opencode/recon-residuals-20260916-161957/` (`evidence/originals/`,
`evidence/*.diff`, `runs/`, `full-vs234/`, `full-sfnet197/`, `work/`)

This closes the two residuals left open by the inertia-apply reconstruction
(§"Inertia-apply reconstruction", §3): `iface_tangref` (target tied to the
pre-`bf74f28` scalar record semantics) and `iface_gen_method_gen` (the
reconstruction had turned the GENERATED interface element into an input element
and dropped the genuine input `element_geometry` tags).

## Result

| corpus | PASS | RUNFAIL | PARSE | total |
|---|---|---|---|---|
| validation-suite 234 (before) | 229 | 5 | 0 | 234 |
| validation-suite 234 (after)  | **231** | **3** | **0** | 234 |
| validation-suite 233 (after, real tests) | **231** | **2** | **0** | 233 |
| sfnet 197 (before / after)    | 187 | 9 | 1 | 197 |

Per-test diff vs `inertia-recon-20260916/full-vs234/vs.tsv`: exactly two
RUNFAIL -> PASS (`iface_tangref`, `iface_gen_method_gen`), **0 regressions**.
The remaining 2 real RUNFAIL are the pre-existing reconstruction artifacts
`elmethod` and `numit` (lost `-follow_material` reference geometry; untouched,
as instructed). The 3rd non-PASS file is the generated mesh dump
`print_mesh_dof.dat` (first line `1 0 0`, not a test input). sfnet 197 is
byte-identical to the stored baseline (187 / 9 / 1).

## 1. `iface_tangref` — honest per-IP value and target recalibration

Since `bf74f28` the history `element_interface_force_tang2` is a length-4
per-integration-point array: `interface.cc` PUTs `f_t2_el_ip`, the elastic trial
`kt2*gamma2_total` of each facing pair (since F4 T6d C52, 2026-10-11, the
clamped force `f_t2_ip`; `iface_tangref` is unchanged) (the `ELEMENT_INTERFACE_FORCE_TANG2`
record, `data_length = 4`). The stored target `10.4025 ±0.2` was set at
`7b414fa` (2026-08-24) from the pre-`bf74f28` **scalar** "side-average,
all-against-all" assembly; the pristine 2026-08-29 `.dbs` reads the scalar
`element_interface_force_tang2 2 10.3532842167`.

Honest-value evidence (HEAD `640855a`, LC_ALL=C, fresh cwd, `runs/fixed-tangref/`):

- HEAD returns `4.99979673431 4.99979673115 4.99979673431 4.99979673115` —
  the four facing pairs agree to 4e-9.
- The model has a single load path. `force_element_volume 0 0 5 0` applies
  `fy = 5` to element 3 only (`force_element_volume_geometry 0 3` is an
  element-number list, `force.cc:66-71`), the unit cube `hex8 9..16`
  (volume 1). Element 2 is the interface (group 10); it is the only bridge
  between element 1 and element 3, and element 1's `x=1` face (nodes 2,4,6,8)
  is fully fixed in y (`bounda_unknown 1`/`2`). The interface therefore carries
  the whole `5` in y.
- The frame is normal `-x`, `t1 = +z` (reference point `1 0.5 5`), `t2 = +y`
  (manual-user `group_interface_tangential_reference_point.md`), so the y-load
  lands entirely in `_tang2`; `_tang` is `-1.75e-14` (target `0 ±0.001`, met).
- The stored per-IP value equals the interface total because the IP weights are
  normalised (`sum w_ip = 1`; uniform `1/4` for the 3D `hex8`,
  `interface.cc:373-374`): `sum w_i * f_i = 4.9998`. Equilibrium with the
  applied `5` confirms `4.9998` is the honest per-IP value.
- The intermediate per-IP values of the bisect (`bf74f28` ~50.8, `4b6548a`
  ~49.5) predate the `6cafafe` interface-measure + honest-Bi-CG changes; they
  are not the converged equilibrium of the current per-pair assembly. Only the
  post-`6cafafe` `4.9998` matches the applied load.

Recalibration: `target_value  2  10.4025 0.2` -> `target_value  2  4.9998 0.2`
(tolerance kept). `rc=0`; the first-IP value `4.99979673431` meets the target.

## 2. `iface_gen_method_gen` — the generated element was mis-classified as input

The reconstruction's drop-list (`CORPUS-RECONSTRUCTION.md` §Method step 2)
classified `element_geometry` as a run-produced record and dropped all of them,
while keeping every `element` record of the dump as input. For this test that
produced a runnable-looking input with the GENERATED interface element 3 as an
input `quad4` and no `element_geometry` at all -> the target `element_geometry 3`
does not exist -> "Error in calculation".

Discriminator from the archived `.dbs`: the **`element_group` record count**.
The archive has three elements but only two `element_group` records:

```
element  1  -quad4 1 2 3 4
element  2  -quad4 5 6 7 8
element  3  -quad4 2 4 5 7        <- generated: NO element_group
element_group  1  0
element_group  2  0
element_geometry  1  20
element_geometry  2  30
element_geometry  3  10           <- generated: method_generate = -element_geometry
```

`control_mesh_generate_interface_method 0 -element_group -element_geometry`
selects by `element_group` and, for the generated element, writes an
`element_geometry` record INSTEAD of an `element_group` one (manual-user
`control_mesh_generate_interface.md` §"Selección y generación por
`element_geometry`"). Element 3 has no `element_group` and carries
`element_geometry 3 10` (= the interface group 10), so it is generated.
`element_geometry 1 20` / `2 30` are genuine **input** tags of the two blocks
(the generator only touches the generated element). The sibling
`iface_gen_method` confirms the pattern: its archive has `element_geometry 1 20
/ 2 30` as input and `element_group 3 10` for its generated interface.

Fix: remove the input `element 3` line and restore the input `element_geometry 1
20` / `2 30` lines. Result (`runs/fixed-gengen/`): `rc=0`; the generator
produces element 3 (`quad4 2 4 5 7`) with `element_geometry 3 10`, and the
target `element_geometry 3 = 10` now reads the **generated** record. The run's
`element_geometry` set is byte-identical to the archive (`1 20 / 2 30 / 3 10`,
`element_group 1 0 / 2 0`), so the test exercises the intended generation
feature, not the input tag.

Family note (NOT changed here — five currently-passing tests would be touched):
the same drop-list defect turned the generated interface element(s) into input
elements in `iface_gen` (element 3), `iface_gen_method` (element 3),
`iface_gen_geom` (element 3) and `iface_gen_quad9` (elements 3 and 4). Their
targets are `node_dof` values, so they pass either way, but their reconstructed
input does not exercise the generation feature. Correcting them is a separate
work unit and is recorded here for follow-up.

## 3. What changed (exact corpus edits)

The corpus directory is gitignored/untracked, so the edits are recorded here
verbatim. Originals are preserved under `evidence/originals/` (sha256 in
`evidence/originals.sha256`; unified diffs in `evidence/*.diff`).

```
iface_tangref.dat:        target_value  2  10.4025 0.2   ->  target_value  2  4.9998 0.2
iface_gen_method_gen.dat: - element  3  -quad4 2 4 5 7
                          + element_geometry  1  20
                          + element_geometry  2  30
```

No other line of either file changed. No source file was modified and no other
RUNFAIL test was touched.

## 4. Reproducibility

- Baseline (`640855a`, before edits): validation-suite 234 `229 / 5 / 0`
  (`inertia-recon-20260916/full-vs234/vs.tsv`); both residuals `rc=1`
  (`runs/base-tangref/`, `runs/base-gengen/`).
- After the edits, same runner and binary: validation-suite 234
  `231 / 3 / 0`; 233 real tests `231 / 2 / 0` (`full-vs234/vs.tsv`).
- sfnet 197: `187 / 9 / 1`, classification byte-identical to the baseline
  (`full-sfnet197/sfnet.tsv`).
- Candidate A/B for `iface_gen_method_gen` (`work/candidates/`):
  `C1` = element 3 removed, no input geometry tags -> `rc=0`, generated
  `element_geometry 3 10`; `C2` = C1 + restored `element_geometry 1 20 / 2 30`
  -> `rc=0`, byte-identical `element_geometry` set to the archive (chosen).

---

# `-follow_material` reference geometry + initial truss force — the last two RUNFAIL (2026-09-16)

Closes `elmethod` and `numit`, the two remaining validation-suite RUNFAIL.
validation-suite 233 real tests: **233 / 0 / 0**; the 234th file
(`print_mesh_dof.dat`) is a generated dump, not a test. sfnet 197 unchanged.

## Result

| corpus | PASS | RUNFAIL | PARSE | total |
|---|---|---|---|---|
| validation-suite 234 (before) | 231 | 3 | 0 | 234 |
| validation-suite 234 (after)  | **233** | **1** | **0** | 234 |
| validation-suite 233 (real tests) | **233** | **0** | **0** | 233 |
| sfnet 197 | 187 | 9 | 1 | 197 |

`elmethod` / `numit` RUNFAIL -> PASS, 0 regressions; the single non-PASS file
is the generated dump `print_mesh_dof.dat`.

## 1. Two independent input records were lost, not one

The previous note (`CORPUS-RECONSTRUCTION.md` §3) said only the
`-follow_material` reference geometry (`node_start_refined`) was lost and that
re-inserting it did not help. That attempt was **wrong in two ways**:

1. it *added* a `node_start_refined` record while keeping the deformed dump
   `node` as the input geometry — so the truss was already at its final shape
   and nothing moved;
2. it missed the second, decisive record: the **initial `element_truss_force`
   self-stress**, dropped by the reconstruction drop-list.

### 1a. The dump `node` is the DEFORMED geometry; `node_start_refined` is the input reference

Direct evidence from a surviving original `.dat`/`.dbs` pair that uses
`-follow_material`, `dpline_method`:

```
dpline_method.dat (input)      : node 3 0 0.98      node 4 1 0.98
dpline_method.dbs (dump)       : node 3 0 0.96      node 4 1 0.96   <- deformed
dpline_method.dbs (dump)       : node_start_refined 3 0 0.98
                                 node_start_refined 4 1 0.98       <- the input
```

`top.cc:269` copies the input `NODE` into `NODE_START_REFINED` at the start of
the run (`if ( max<0 ) db_copy( NODE, NODE_START_REFINED, ... )`), and under
`-FOLLOW_MATERIAL` the mesh then moves (`truss.cc:107`, `locate.cc:47`), so the
final `node` record is deformed while `node_start_refined` keeps the reference.
The reconstructed `elmethod`/`numit` used the **deformed** `node`
(`0.877` / `1.631`) as input; the correct input is the reference
(`0` / `1` / `2`), i.e. the values the dump stores in `node_start_refined`.

### 1b. `element_truss_force` is a genuine INPUT record and was dropped

- It is a registered data item (`database.cc:2918`: `ELEMENT`, `DOUBLE_PRECISION`,
  `data_length = 1`), so it is parseable in the data part.
- It is in the reconstruction drop-list (`reconstruct.py:46`), so the
  reconstruction removed it from every test.
- The dump stores the **final** (relaxed) value, `element_truss_force ≈ 0`, so
  the input self-stress is not visible in the `.dbs` — it must be derived.
- `ProjectDocs/manual-developer/control_print_element_method.md` §Detalles
  documents the mechanism: *"la fuerza de un `element_truss_force` inicial es un
  self-stress que SE RELAXA con la deformación (`new_force = old_force + k*dL`;
  en equilibrio final queda 0)"*.

### 1c. Derivation of the lost initial forces (evidence-led, unique)

`truss.cc:181-188`: `new_truss_force = old_truss_force + truss_stiffness *
incremental_length`. With `group_truss_memory 0 -updated_without_rotation` the
`else` branch runs (`memory==-UPDATED` is a different enum), so
`truss_stiffness = group_truss_young * group_truss_area / initial_length`
(`truss.cc:160`). Here `E = A = 1`, `initial_length = 1` (the reference
geometry), so `k = 1`. At the archived equilibrium `new_truss_force = 0`, hence

```
F_init = k * (initial_length - final_length)
```

The archived geometry (`node` final, `node_start_refined` reference) gives the
final lengths, so the initial forces are fully determined:

| test | element | L_ref | L_final | F_init |
|---|---|---|---|---|
| elmethod | 1 (nodes 1-2) | 1 | 0.877 | **0.123** |
| elmethod | 2 (nodes 2-3) | 1 | 1.631 - 0.877 = 0.754 | **0.246** |
| numit | 1 (nodes 1-2) | 1 | 0.877 | **0.123** |

A single `group_truss_initial_force` cannot produce this: the same group value
would give `disp(node2) = -F0` and `disp(node3) = -2*F0`, i.e. a ratio of 2,
whereas the archive has `-0.123` / `-0.369` (ratio 3). The forces are
per-element `element_truss_force` records (both elements are in group 0: the
archive has no `element_group` record).

## 2. Verification — the reconstruction reproduces the archive exactly

`elmethod` (fixed input, current binary, fresh cwd, `LC_ALL=C`):

```
node  1  0          node_dof  2  -2.22044604925e-15 -0.123 0 0 0 0 0 0
node  2  0.877      node_dof  3  -4.4408920985e-16 -0.369 0 0 0 0 0 0
node  3  1.631      element_truss_force  1  -1.11022302463e-16
node_start_refined  1 0 / 2 1 / 3 2       2  1.11022302463e-16
```

This is **byte-identical** to the archived `.dbs` (including the solver noise
`-2.22e-15` / `-4.44e-16`), and the print files match too:
`element_truss_force_0.20` = `0.4385 1.11022302463e-16 / 1.254 -1.11e-16`,
`element_truss_force_0.21` = `0 -1.11e-16 / 0.877 -1.11e-16 / 0.877 1.11e-16 /
1.631 1.11e-16`. `numit` likewise gives `node_dof 2 = 0 -0.123 ...` exactly as
the archive, with exactly 16 `control_print_number_iterations` monitor lines.

The primary features are unaffected: `elmethod` still discriminates
`-middle` (2 lines) vs `-node` (4 lines) with `2*mid == x_node2`; `numit` still
prints 16 monitor lines. The target is met (`rc=0`) for both.

## 3. What changed (exact corpus edits)

The corpus directory is gitignored/untracked, so the edits are recorded here
verbatim. Originals are preserved under `evidence/originals/` (sha256 in
`evidence/originals.sha256`; unified diffs in `evidence/*.diff`).

```
elmethod.dat: - node  2  0.877
              - node  3  1.631
              + node  2  1
              + node  3  2
              + element_truss_force  1  0.123
              + element_truss_force  2  0.246

numit.dat:    - node  2  0.877
              + node  2  1
              + element_truss_force  1  0.123
```

The `element_truss_force` lines are inserted right after the `element` records
(their canonical position in the dump order). No other line changed, no source
file was modified, and no other test was touched.

## 4. Rejected hypotheses (recorded, not applied)

- **Add `node_start_refined` only** (the previous attempt): keeps the deformed
  `node`, so `initial_coord = coord - old_dof = coord` and the truss is already
  at rest -> target still `0` (`runs/C3-*`).
- **Seed `node_dof -velix` with the target** (`-0.123` / `-0.369` as an initial
  condition): makes `rc=0` but is circular (it seeds the answer) and leaves
  `node_start_refined` = the deformed geometry, contradicting the archive
  (`runs/C4-*`).
- **`group_truss_initial_force`**: cannot produce the observed 1:3 displacement
  ratio (see §1c).

## 5. Reproducibility

- Binary: `build/tochnog` rebuilt with `scripts/build_safe.sh --clean`
  (HEAD `92518b8`, after the `top.cc` stress-pass hoist `3957809`); suite
  build_safe 16/16 runs OK.
- Env: `LC_ALL=C`, `LD_LIBRARY_PATH=<repo>/external-downloads/numlib-runtime`,
  fresh cwd per test, `timeout 45` (runner `work/run-corpus.sh`).
- validation-suite 234 after: `233 / 1 / 0` (`full-vs234/vs.tsv`); the only
  RUNFAIL is `print_mesh_dof` (generated dump).
- sfnet 197: `187 / 9 / 1`, classification byte-identical to the baseline
  (`full-sfnet197/sfnet.tsv`).
- Candidates: `runs/baseline-*` (deformed node, target `0`), `runs/fixed-*`
  (reference node only, target `0`), `runs/C3-*`, `runs/C4-*`, `runs/C5-*`
  (reference node + initial force = the archive).
