# solver_matrix_symmetric

## Where implemented

- `so_bicg.cc`, `solve_iterative_bicg()` (per-solve solver selection,
  around lines 525-555): with `solver_matrix_symmetric -yes` the GNU
  previously forced plain CG even when the per-solve symmetry
  measurement of the assembled system failed, which diverges on coupled
  velocity-pressure systems. Now the honest Bi-CG path runs on the
  as-assembled matrix when the measurement disagrees with the user's
  assertion (warning printed), and plain CG when the system measures
  symmetric.

## Implementation details

The per-solve measurement (`solve_iterative_bicg_symmetric`) checks the
assembled element matrices; the choice between CG and Bi-CG is
`solve_iterative_bicg_use_cg = measured_symmetric` for both values of
the switch, so `-yes` only adds the warning. The Professional (manual
6.1052) symmetrizes the matrix "if needed" and then runs a symmetric
solver; the GNU does not symmetrize: Bi-CG on the true system is used
instead (the anti-symmetric part is honoured rather than averaged
away). Pure-mechanical and bounda_alternate-staggered problems measure
symmetric and keep the CG path (identical results).

No external library; no hardcoded parameters beyond the symmetry
measurement of `solve_iterative_bicg_symmetric`.

## Verification

- `large2` (3D coupled consolidation brick under a surface load, pres
  drainage at the far face): previously diverged (CG on a
  non-symmetric system, residual grew over 49000 iterations); now runs
  to completion with one transient Bi-CG non-convergence recovered by
  the direct-LU retry. rc=0 (no targets; solver/memory test).
- `large3` (same model + bounda_alternate + `control_timestep_
  iterations`): the alternated systems measure symmetric; unchanged,
  rc=0.
- K24 guard tests (F4 T6b, 2026-10-11), `validation-suite/test-2026`,
  wired in `scripts/build_safe.sh` (`run_ok`), targets = Professional
  02-08-2026 with `-no` at node 1 (0, -0.5, -0.5):
  - `k24_consol_3d1` (`-yes`) / `k24_consol_3d2` (`-no`): transient,
    k 0.1, dt 0.1 to t 2.2: pres -0.7085296 (tol 7e-4), disx 0.05420111
    (tol 5.4e-5); GNU -0.708335 / 0.0542098 for both switches.
  - `k24_consol_3d_undrained1` (`-yes`) / `_undrained2` (`-no`): k 1e-8,
    one step: pres -0.847677 (tol 8.5e-4), disx 0.01748967 (tol 1e-4);
    GNU -0.8476521 / 0.0175445 for both switches.
  The Professional `-yes` values (-0.678332 / disx 0.05644; -0.803176 /
  disx 0.02077) fail these targets.

## Differences with the Professional version

- K24 (DIVERGENCES-PROFESSIONAL.md, class A, project decision
  2026-10-11; also Pro-suspect): the Professional symmetrizes the
  non-symmetric Biot u-p matrix under `-yes` and its pore pressure is
  non-physical (undrained drift with no flow, wrong sign for small
  capacities, 4.5 % off on `large2`); the GNU keeps the as-assembled
  matrix and equals the Professional run with `-no`. Evidence in the
  user page and in the F4 T6b entry of `odd/tasks/f4-investigations.md`.

## Pending

- Symmetrization of the assembled matrix (Professional parity) is not
  implemented, deliberately (K24).
- The stale `so_bicg.cc` comments next to the selection (CG under
  `-yes`, manual 6.1053) were removed/corrected after F4 T6b.
