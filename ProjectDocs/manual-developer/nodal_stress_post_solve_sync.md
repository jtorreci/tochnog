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

### `top.cc` `top()` — post-solve pass

Immediately after `solve(options_solver)` and before
`parallel_new_dof_diagonal()`, when `stress_recovery_absolute && materi_stress
&& stres_indx>=0 && !find_local_softvar`:

1. Save the full `NODE_LHSIDE`/`NODE_RHSIDE` vectors of every active node.
2. Zero the six stress slots (`ip0 = stres_indx/nder`, `ip0..ip0+5`).
3. Set `stress_recovery_post_pass = 1`, call `element_loop()`, clear it.
4. Restore every **non-stress** entry from the saved copy.

Step 4 is required because structural elements (`truss`, `beam`, `spring`,
`contactspring`, interface) bypass `materi()` via `goto
add_to_system_vectors`; without the restore their first-pass contributions
would be doubled (measured: `truss2`, `truss3`, `conspr4` regressed before the
restore was added). Continuum elements are protected by the `materi()` gate
below, so the restore is a no-op for them.

`parallel_new_dof_diagonal()` then runs once per iteration with a mix of:
fresh stress entries (from the solved velocity) and the pre-solve non-stress
entries. It also recomputes the `materi_stress_pressure_history` dof from the
refreshed stress slots.

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

## Verification

- `mpc5.dat`: `sigyy = 0.5` at `k=1, 2, 3, 5` (target `0.5 1e-3`, rc=0).
- Family `delete1..4`, `mpc1..9`: all rc=0; `delete2` node 2 `sigxx = 0.8`,
  `delete3` `sigyy = 0.6`, byte-identical to the baseline.
- sfnet 197: **187 PASS / 9 RUNFAIL / 1 PARSE**, per-test classification diff
  vs baseline = **0**.
- validation-suite 233: **201 PASS / 32 RUNFAIL / 0 PARSE**, diff = **0**.
- Internal suite: **16/16** runs OK.

## Pending refactorings

- The pass currently runs on every equilibrium iteration. Only the **last**
  iteration's report is observed, so it could be hoisted after the loop
  (updating the stress dofs manually or via a stress-only variant of
  `parallel_new_dof_diagonal`), halving the added cost. Not done because the
  per-iteration form also keeps any in-loop print in sync.
- The save/restore of the whole nodal RHS/LHS vectors is memory-heavy; it
  could be reduced to the non-stress slots only.
- A cleaner long-term fix is the monolithic saddle-point mixed solve
  (`DIAG-SOLVE-MIXTO.md` §8 C), which would make the stress principal and
  remove the lag structurally; that remains refuted as structurally
  inconsistent for shared nodes in the current tie scheme.
