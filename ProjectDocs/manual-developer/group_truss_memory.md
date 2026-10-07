# group_truss_memory — developer notes

## Where

`truss.cc`, `truss()`: the "direction + deformations + stiffness" block
(memory branches) and the element matrix loop.

## K17 fix (2026-10-06)

- **`-total` branch** (new, `memory==-TOTAL && truss_reference!=-YES`):
  direction from the CURRENT coordinates (`new_coord`), `new_length` /
  `old_length` the current lengths, `incremental_length = new - old`,
  `truss_stiffness = EA / initial_length`. Accumulated with the stored
  `ELEMENT_TRUSS_FORCE` this is `N = N0 + EA (L - L0)/L0`, the force the
  Professional 02-08-2026 reports on `truss14` (N = 1.414214 at L = 2.414).
  Before, `-total` fell into the linear `else` branch (initial axis), so a
  transverse load met no stiffness (bicg breakdown, disy 1e10).
- **Geometric stiffness** (`geometric_nonlinear` flag, set by the
  `-updated` and `-total` branches): `tmp += fac*dt*(N_mix/L)*(delta_ij -
  n_i n_j)` in the matrix loop, only when `N_mix > 0`. `N_mix` is the
  force the RHS uses (materi_dynamic blend), so the term is the exact
  derivative of `fac * n * N_mix` with respect to the transverse velocity.
  RHS untouched: converged states do not change (validation-suite /
  test-2026 dumps identical except `bmom_truss`, 3.6e-7; sfnet `truss2/3`
  2.9e-5 and `truss7` 3e-5 relative, iteration-path differences within
  the automatic-iteration tolerance; no exit code changes).
- **Hardcoded**: tension only (`N > 0`). A compressive N would lower the
  transverse pivot of a free node below zero and the bicg would break
  down; the consistent tangent in compression is left out on purpose.
- **Not changed**: the `truss_reference` path (displacement-primary
  dynamics) keeps the linear reference geometry; the zero-force first
  iteration of a straight truss is still singular (bicg "broke down"
  message, run continues; BACKLOG C12). Pro solves it with its default
  (direct) solver; Pro with `-matrix_iterative_bicg` also stops on truss2.

## Verification

- `truss_total_rotation1` (= Pro truss14 with targets): RED rc 1 (disy
  1e10), GREEN disx 0.707107 / disy 1.707107 / N 1.414214 (Pro exact).
- `truss_updated_transverse1` (= Pro truss2 plus targets): RED rc 1
  (reaction velx -0.99997), GREEN reactions 0.1/0.1, veliy 0.823578 (Pro
  0.823577).
