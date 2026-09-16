# nodal_stress_post_solve_sync

Status: **IMPLEMENTED 2026-09-15**. Closes the last documented `mpc5`
residual: the reported nodal stress no longer trails the solved velocity by
one iteration. `mpc5.dat` now reports `sigyy = 0.5` at `k=1` (and still `0.5`
at `k=2`, matching the Professional binary 02-08-2026).

## Symptom

`mpc5.dat` (quad9 partially deleted via
`control_mesh_delete_geometry_factor`, with a `control_reset_dof -sigyy` of 1
and a non-conforming `mpc_linear_quadratic` tie) expects
`post_point_dof 10 -sigyy = 0.5` (target `0.5 1e-3`).

| k | GNU before | GNU after | Professional |
|---|---|---|---|
| 1 | 1.0 | **0.5** | 0.5 |
| 2 | 0.5 | **0.5** | 0.5 |
| 3 | 0.5 | 0.5 | 0.5 |

The `k=2` value already matched the Professional after the absolute-weighted
nodal-stress recovery (`nodal_stress_recovery`, commit `0b23957`); the `k=1`
report was the remaining lag.

## Root cause

The `materi_stress` dofs are **non-principal** (`input.cc`: no
`dof_principal` assignment). They are therefore not part of the solved
unknown vector; `dof.cc parallel_new_dof_diagonal` updates them **after** the
solve as

```
node_dof_new[stres] = node_rhside[stres] / node_lhside[stres]
```

using the `element_rhside`/`element_lhside` assembled in `element_loop()`
**before** the solve. `materi.cc` assembles the recovery numerator as
`volume * weight * new_sig`, where `new_sig` is the constitutive stress
computed from the **pre-solve** `new_unknowns` (the previous iterate's
velocity). The reported stress therefore lags the velocity by one iteration.

## Candidate (b) verification (measured)

The scheduled note (`mpc5_staggered_stress_lag.md`) listed two candidates:
(a) a monolithic saddle-point solve and (b) a post-solve constitutive
re-evaluation. A previous worker argued that (b) is "algebraically identical
to running iteration k+1".

Measured verdict: the **value** produced by candidate (b) at iteration `k`
equals the value the unmodified scheme reports at iteration `k+1`
(`proj(sigma_constit(v_k))`), but the **mechanism is distinct**: candidate (b)
adds only a constitutive evaluation (no extra solve), so it syncs the `k=1`
report without running a second momentum solve. The claim "identical to
iteration k+1" is true of the value sequence, not of the cost, and it does
**not** prevent the `k=1` sync. The recovery fix (`0b23957`) made the `k=2`
value correct (0.5), so the one-iteration-ahead report at `k=1` is now exactly
the Professional's `0.5`.

The concept was first confirmed with a runtime probe (a full second
`element_loop()` after `solve()` gated by `TOCHNOG_STRESS_POST`): `k=1` went
from `1.0` to `0.5`. The production implementation below is the same idea
restricted to the stress recovery.

## Implementation

### Global flag

- `long int stress_recovery_post_pass = 0;` — defined in `polynom.cc`
  (next to `stress_recovery_absolute`), declared `extern` in `tochnog.h` and
  `tochnog-mod.h` (the two headers must stay in sync).

### `top.cc` `top()` — post-solve pass (hoisted after the equilibrium loop)

**Since 2026-09-16** the pass runs once per time step, immediately **after**
the equilibrium `for ( iteration ... )` loop and before
`db_version_copy(VERSION_NEW, VERSION_NORMAL)` (i.e. before the step close),
when `nuknwn>0 && stress_recovery_absolute && materi_stress && stres_indx>=0
&& !find_local_softvar`:

1. Save the full `NODE_LHSIDE`/`NODE_RHSIDE` vectors of every active node.
2. Zero the six stress slots (`ip0 = stres_indx/nder`, `ip0..ip0+5`).
3. Set `stress_recovery_post_pass = 1`, call `element_loop()`, clear it.
4. Restore every **non-stress** entry from the saved copy.
5. Assign the stress dofs from the refreshed recovery: for each active node
   `node_dof_new[stres_indx+k*nder] = node_rhside[ip0+k]/node_lhside[ip0+k]`
   (`k=0..5`, when the slot is not bounded and `lhside!=0`), then refresh the
   `materi_stress_pressure_history` (`sph`) dof and, when `derivatives` is
   active, the `tsig*` stress-rate slots, exactly as
   `dof.cc parallel_new_dof_diagonal` does. Only the stress slots are touched:
   re-running the whole function would double the `+=` accumulation of the
   other non-principal dofs.

Step 4 is required because structural elements (`truss`, `beam`, `spring`,
`contactspring`, interface) bypass `materi()` via `goto
add_to_system_vectors`; without the restore their first-pass contributions
would be doubled (measured: `truss2`, `truss3`, `conspr4` regressed before the
restore was added). Continuum elements are protected by the `materi()` gate
below, so the restore is a no-op for them.

The in-loop `parallel_new_dof_diagonal()` still runs every iteration; during
the iterations it assigns the stress from the pre-solve assembly, and the
hoisted block overwrites the stress dofs with the refreshed final recovery.

The pre-hoist form (kept here as the design record) ran the same 1-4 block
between `solve()` and `parallel_new_dof_diagonal()` on **every** iteration.
The hoist is behaviourally equivalent for the required targets (mpc5/delete2/
delete3, sfnet 197, validation-suite 233); two non-target derived-dof changes
are documented in the manual-user note and in `IMPROVEMENT-ANALYSIS.md` item
(a): the `xsigyy` x-derivative of `contact_heatgen` (the x/y stress rates come
from the previous iteration's stress in both forms) and the `mpc8` `sigxx`
field (`-total_linear` + `materi_displacement`: the hoisted pass sees the
updated displacement, so the reported stress is consistent with the reported
displacement).

### `materi.cc` `materi()` — stress-only assembly

When `stress_recovery_post_pass` is set:

- the velocity/momentum block is skipped (`if (!stress_recovery_post_pass)`);
- the `h>EPS_H` block before the stress recovery is skipped;
- after the stress-recovery block a `if (stress_recovery_post_pass)
  continue;` skips every other dof block.

Only the stress recovery (`volume*weight*new_sig` numerator and
`volume*weight` denominator) is assembled.

### `elem.cc` `elem()` / `element_loop()` — stress-only pass

When `stress_recovery_post_pass` is set:

- `general()` and the non-`MATERI` type calls (`condif`, `groundflow`,
  `maxwell`, `wave`) are skipped; only `materi()` runs;
- the post-`materi` bookkeeping blocks (derivatives, strain energy,
  `options_element_dof`, `element_dof_tmp`) are skipped;
- the post-integration-point block (tendon normalization, derivatives
  normalization, `area()`, residue, `ELEMENT_VOLUME`/`ELEMENT_MASS`) is
  skipped;
- the sparse matrix store is skipped (the matrix is not needed after the
  solve);
- the deleted-element block does **not** re-apply the frozen released load
  (`element_rhside_delete` is already stored); it only zeroes the stress slots
  of both `element_rhside` and `element_lhside`, so a deleted element
  contributes neither numerator nor denominator to the recovery;
- in `element_loop()`'s merge, the `derivatives` write to `NODE_DOF` is
  skipped.

The final merge into `NODE_LHSIDE`/`NODE_RHSIDE` runs normally (it only sees
the stress entries), which is what makes `parallel_new_dof_diagonal` pick up
the refreshed recovery.

## External features / libraries

None. Pure C++ inside the existing element/assembly pipeline.

## Hardcoded parameters

- The six stress components: `ip0..ip0+5` (the compacted symmetric layout of
  `stress_indx()`: `sigxx, sigxy, sigxz, sigyy, sigyz, sigzz`).
- The pass is gated on `stress_recovery_absolute` (absolute weighted average),
  so convected stress dofs (`materi_strain_total`) keep the historic
  incremental update and are not re-evaluated post-solve.

## Cost

The element pass runs twice per equilibrium iteration for models with
`materi_stress` + the absolute recovery (the second pass assembles only the
stress block). This is the price of reporting in sync; the corpus shows no
timeout or classification change.

**HOISTED 2026-09-16** (commit `dev/stress-pass-hoist`): the pass no longer
runs inside the equilibrium loop. Only the last iteration's report is
observed, so the stress-only `element_loop()` and the nodal stress update run
**once per time step**, after the loop (see `top.cc` below). Measured on
representative `materi_stress` + absolute-recovery models (same machine,
`LC_ALL=C`, fresh cwds, `build_safe.sh`):

| test | elements | equilib. iters | base min/median (ms) | hoisted min/median (ms) | speedup |
|---|---|---|---|---|---|
| mpc5 | 2 | 2 | 17 / 19 | 17 / 18 | startup-dominated |
| delete2 | 2 | 4 | 16 / 18 | 16 / 17 | startup-dominated |
| delete3 | 2 | 25 | 35 / 38 | 31 / 33 | ~1.14x |
| dynamic8 | 16 | 4000 | 11227 / 12180 | 10009 / 10028 | 1.12-1.21x |
| large3 | 784 | 220 | 56972 / 59024 | 46224 / 47213 | 1.23-1.25x |

`large3` is 22 steps x 10 iterations: the stress-only passes drop from 220 to
22. The added per-iteration assembly cost is therefore removed for all but the
last iteration of each step.

## Verification

- `mpc5.dat`: `sigyy = 0.5` at `k=1, 2, 3, 5` (target `0.5 1e-3`, rc=0).
- Family `delete1..4`, `mpc1..9`: all rc=0; `delete2` node 2 `sigxx = 0.8`,
  `delete3` `sigyy = 0.6`, byte-identical to the baseline.
- sfnet 197: **187 PASS / 9 RUNFAIL / 1 PARSE**, per-test classification diff
  vs baseline = **0**.
- validation-suite 233: **201 PASS / 32 RUNFAIL / 0 PARSE**, diff = **0**.
- Internal suite: **16/16** runs OK.

**Hoist verification (2026-09-16, base `640855a` vs `dev/stress-pass-hoist`)**:

- mpc5 `k=1,2,3,5` -> `sigyy = 0.5` rc=0 (base and hoisted identical);
  delete2 node 2 `sigxx = 0.8`, delete3 `sigyy = 0.6`; the `node_dof` +
  `post_point_dof` fields byte-identical.
- sfnet 197: **187 PASS / 9 RUNFAIL / 1 PARSE** in both, per-test
  classification diff = **0**; `node_dof`/`post_point_dof` value diff = **0**.
- validation-suite 233 (excl. the generated `print_mesh_dof.dat`):
  **231 PASS / 2 RUNFAIL** (`elmethod`, `numit`) in both, classification diff
  = **0**; one non-target value change (`contact_heatgen` `xsigyy`).
- Internal suite: **16/16** runs OK for both (full-mode file-check FALLO list
  identical).

## Pending refactorings

- ~~The pass currently runs on every equilibrium iteration.~~ **DONE
  2026-09-16** (hoisted after the loop; see `top.cc` above). Residual: the
  hoisted pass does not refresh the `derivatives` x/y stress rates (they are
  computed by the normal `element_loop` from the previous iteration's stress
  in both forms) nor the `-total_linear` stress of `materi_displacement`
  models from the pre-displacement-update state; both are non-target outputs.
- The save/restore of the whole nodal RHS/LHS vectors is memory-heavy; it
  could be reduced to the non-stress slots only.
- A cleaner long-term fix is the monolithic saddle-point mixed solve
  (`DIAG-SOLVE-MIXTO.md` §8 C), which would make the stress principal and
  remove the lag structurally; that remains refuted as structurally
  inconsistent for shared nodes in the current tie scheme.
