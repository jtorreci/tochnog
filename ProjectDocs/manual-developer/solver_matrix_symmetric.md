# solver_matrix_symmetric

## Where implemented

- `so_bicg.cc` (per-solve solver selection): with
  `solver_matrix_symmetric -yes` the GNU previously forced plain CG
  even when the per-solve symmetry measurement of the assembled system
  failed, which diverges on coupled velocity-pressure systems. Now the
  honest Bi-CG path runs on the as-assembled matrix when the
  measurement disagrees with the user's assertion (warning printed),
  and plain CG when the system measures symmetric.

## Implementation details

The per-solve measurement (`solve_iterative_bicg_symmetric`) checks the
assembled element matrices; the choice between CG and Bi-CG is
`solve_iterative_bicg_use_cg`. The Professional (manual 6.1053)
symmetrizes the matrix "if needed" and then runs a symmetric solver;
the GNU does not re-symmetrize: Bi-CG on the true system is used
instead (the anti-symmetric part is honoured rather than averaged
away). Pure-mechanical and bounda_alternate-staggered problems measure
symmetric and keep the CG path (identical results).

## Verification

- `large2` (3D coupled consolidation brick under a surface load, pres
  drainage at the far face): previously diverged (CG on a
  non-symmetric system, residual grew over 49000 iterations); now runs
  to completion in ~19 s with one transient Bi-CG non-convergence
  recovered by the direct-LU retry. rc=0 (no targets; solver/memory
  test).
- `large3` (same model + bounda_alternate + `control_timestep_
  iterations`): the alternated systems measure symmetric; unchanged,
  rc=0 in ~34 s.
- Regression surface: only the `large*` corpus tests use
  `solver_matrix_symmetric`; the GNU suite has no other user.

## Pending

- Actual symmetrization of the assembled matrix (Professional parity)
  is not implemented; the honest Bi-CG fallback converges on the same
  problems without touching the solver core.
