# mpc5 — staggered u-sigma stress lag (scheduled work unit)

Status: **PARTIALLY CLOSED 2026-09-15** — the RECOVERY FORM is aligned
(`nodal_stress_recovery.md`: the nodal stress is now the Professional's
absolute weighted average; the full `mpc5` `node_dof` field matches the Pro
to 5e-17 and the target `post_point_dof 10 -sigyy = 0.5` passes at the
test's own 2 iterations). The one-iteration REPORTING lag remains: `mpc5`
still gives `1.0` at `k=1` and `0.5` at `k=2` while the Professional gives
`0.5` at `k=1`, because the stress dofs are non-principal and are updated
after the solve from the PRE-solve constitutive stress. Closing it needs
candidate (b) below (post-solve constitutive re-evaluation); the recovery
weight is no longer the blocker. The rest of this note is the original
design record.

Status (original): **SCHEDULED / BLOCKED** (not implemented). This note
records the precise blocker, the candidate fixes, the estimated scope and
the acceptance criteria so the convergence line can be planned. It is a
design note, not a description of landed behaviour.

## Symptom

`mpc5.dat` (quad9 partially deleted via `control_mesh_delete_geometry_factor`,
with a `control_reset_dof -sigyy` of 1 and a non-conforming `mpc_linear_quadratic`
tie) expects `post_point_dof 10 -sigyy = 0.5` (target `0.5 1e-3`).

With the current GNU mixed u-sigma scheme the run is `rc=1`, `sigyy = 0.7` at
the test's own (default) iteration count. Sweeping
`control_timestep_iterations k`:

| k | GNU | Professional |
|---|---|---|
| 1 | 1.0 | 0.5 |
| 2 | 0.7 | 0.5 |
| 3 | 0.59 | 0.5 |
| 5 | 0.5211 | 0.5 |
| 8 | 0.5026063 | 0.5 |
| 20 | 0.5000006 | 0.5 |

The Professional is one-pass (`0.5` at `k=1`). The GNU relaxes geometrically
with an **asymptotic rate of exactly 0.5**, i.e. the spectral radius of the
staggered u-sigma fixed-point operator for this configuration is 0.5.

## Root cause (measured)

The stress dofs are **non-principal**: `input.cc:525-529` sets no
`dof_principal` for `materi_stress`. A non-principal dof is not part of the
solved unknown vector; it is updated **after** the solve by
`dof.cc parallel_new_dof_diagonal` as
`node_dof_new[stres] += node_rhside / node_lhside`, using the
`element_rhside` that was assembled **before** the solve.

`materi.cc:1042-1044` assembles that RHS as `(new_sig - old_sig)/dt`, with
`new_sig` computed from the **pre-solve** `new_unknowns` (the previous
iterate's velocity). Consequently:

1. the stress lags the velocity by one outer iteration, and
2. the momentum equation in turn uses the previous iterate's stress.

The resulting iteration operator has spectral radius 0.5 for mpc5. This is a
property of the staggered u-sigma fixed point, **not** of the
`element_delete_factor` matrix scaling and **not** of the tie weights: zeroing
the deleted element matrix (released RHS kept) still gives rate 0.5, only a
different fixed point (`1.0, 0.625, 0.5625, 0.515625, 0.50048828125` for
`k = 1,2,3,5,10`).

Note: an earlier "reset-only one-pass" observation was the trivial fixed point
`sigma* = 1`. A minimal mixed probe (1 quad4, `materi_stress`, nodal force)
gives `sigyy = 0.0` at `k=1` and `1.73205080757` (converged) at `k=2`, so the
GNU mixed u-sigma scheme is a 2-pass staggered scheme for any load (GNU default
`iteration_min = iteration_max = 2`, `top.cc:350-353`).

## Candidate fixes

### (a) Monolithic mixed solve — DIAG-SOLVE-MIXTO §8 C

Assemble the true saddle-point system `[0 B^T; B -C^-1]` (remove `dt*K` from
the momentum block when `materi_stress` is active), make the stress dofs
principal and add the coupling (`-D B dt` in the stress row, `B^T` in the
momentum row), and solve with MINRES (Benzi-Golub-Liesen) or SuperLU with
pivoting (`-matrix_superlu` is already linked and working).

- Scope: **large**, corpus-wide. Refactor of `elem.cc` / `so.cc`, new solver
  path, new validation for the whole suite.
- Benefit: correct physics and a well-posed problem; closes the whole
  staggered-family residual, not just mpc5.

### (b) Post-solve constitutive re-evaluation

After the solve, recompute the stress from the solved velocity (a second
element pass) so the momentum equation of the next outer iteration sees a
stress consistent with the current velocity. This is effectively what
iteration 2 does today; making it an explicit, always-on pass would collapse
the fixed point.

- Scope: medium; touches the element pass ordering and the RHS assembly.
- Risk: the second pass is not free; it must not double-count the
  `element_rhside` or the deleted-element frozen RHS.

## Estimated scope

Weeks (candidate (a) is the larger, corpus-wide refactor; candidate (b) is
narrower but needs careful A/B across the staggered family).

## Acceptance criteria

- `mpc5.dat` `rc=0` with `post_point_dof 10 -sigyy = 0.5` (target `0.5 +/- 1e-3`)
  at the test's own default iteration count (the input sets no
  `control_timestep_iterations`), i.e. the fixed point must be one-pass like
  the Professional (`0.5` at `k=1`), not a geometric relaxation.
- Corpus no regression: sfnet 197 classification unchanged from the current
  baseline (160 PASS / 37 RUNFAIL / 0 PARSE) and internal suite 16/16.

## Related

- `ProjectDocs/DIAG-SOLVE-MIXTO.md` §8 C (proposed monolithic mixed solve),
  §15-16 (the value-constrained vs energy-consistent tie history).
- `ProjectDocs/manual-developer/delete_factor_element_matrix_scaling.md` — the
  delete-factor scaling that closed delete3; mpc5 is a distinct mechanism.
- `ProjectDocs/SEGUIMIENTO-CONVERGENCIA.md` — the scheduling row and the mpc5
  evidence references.
