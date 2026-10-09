# Solver memory — developer notes (F4 T7)

## Where it is implemented

- `so.cc`, function `solve()`: workspace allocation (after the eigen
  checks), the band/sparse fill loops and the Bi-CG failure branch.
- `so.cc`, static `solve_band_retry()` and
  `solve_band_retry_entries()`: the lazy band LU retry (B).

## Root cause (A): unused dense matrix

Commit `5b7bb8b` (2026-08-30) added a dense `n x n` matrix `mat_dense`
(plus `ipiv_dense`) to every solve, allocated with `get_new_dbl` and
zero-filled with `array_set`, and filled it in the element and node loops
"for the Bi-CG retry". Its only reader was the `else` branch of the retry
(`dgesv_` when `mat` is NULL), but the band matrix `mat` was always
allocated when `n > 0`, so that branch never ran. Cost per solve: 8 n²
bytes and an O(n²) zero-fill, independent of the model.

Fix: `mat_dense`, `ipiv_dense`, their fills, the dead `dgesv_` branch and
their frees are removed. The retry keeps the band LU (`dgbsv_`). The
`dgesv_` declaration stays (unused) for a future dense route.

Verification: build_safe 154 OK / 0 FALLO; corpus A/B on omen (old
`225d536` vs A) — see the SEGUIMIENTO row of the commit.

## Root cause (B): band workspace on every Bi-CG solve

Since `5b7bb8b` the band LU workspace (`mat`, `(2*kl+ku+1)*n` doubles with
`kl = ku = band`) was also allocated, zero-filled and filled on every
Bi-CG solve, for the direct-LU retry after a Bi-CG failure. On the Bi-CG
path the band came from the raw node numbering (`(1+diff)*nprinc` over
`NODE_NODE`, no reordering), so a mesh with distant node numbers (3D
meshes, `control_mesh_macro`, refined meshes) needed `(3*band+1)*n`
doubles: several GB for `force11` / `large1`. The retry runs in a small
minority of solves.

Fix: `so.cc` `solve_band_retry()` (static, above `solve()`), called only
after `solve_iterative_bicg()` fails:

1. Orders the active nodes with a Cuthill-McKee front (breadth-first over
   `NODE_NODE`, a new front at the first unordered active node) and maps
   every local dof to a band row (`perm`). The Bi-CG ordering
   (`solve_global_local`) is not changed.
2. `solve_band_retry_entries(0, ...)` visits the element matrix entries
   (`ELEMENT_MATRIX_VALUES`, or `GROUP_MATRIX_VALUES` with
   `options_matrix_group`), with the mpc tie expansion of the slave
   elimination, and measures the EXACT band of the node-number ordering
   and of the Cuthill-McKee ordering; the smaller one is used (the
   node-number ordering on a tie, which reproduces the old retry bit for
   bit: the entries outside the exact band are zero).
3. Allocates `mat` with that band, fills the entries (mode 1), adds
   `NODE_LHSIDE` and the `EPS_MAT` diagonal floor exactly as the old fill,
   permutes the original right-hand side (`solve_b_temp`), calls `dgbsv_`
   and un-permutes the correction into `solve_b`.

The element matrices are still in the database at that point (the Bi-CG
only reads them), so the lazy fill reproduces the system of the
iterative solve. The band storage in `solve()` is now allocated only for
the band routes (`-matrix_lapack`, eigen), and the `NODE_NODE` band loop
of the Bi-CG branch is gone. C12/C13 (`bicg_stop`) are unchanged: they
act inside `solve_iterative_bicg()`, before the retry decision.

Results: when the Cuthill-McKee band is smaller the LU pivots in another
order, so retried solves change at round-off level (corpus A/B: `force9`
4e-7, `pressed_rubber_cylinder` 1e-7, `force5` 1e-9 of the column
maximum; `examp2` -0.634519741323 unchanged).

Verification: examp2 (retry on every solve) bit-identical; build_safe (A+B+C
build) 158 OK / 0 FALLO; corpus A/B on omen (A -> A+B): only round-off in the retried
solves listed above, and `force11` now starts solving (bad_alloc ->
Bi-CG over the 120 s corpus slot).

## External features / libraries

- LAPACK `dgbsv_` (band LU).

## Hardcoded parameters

- `EPS_MAT` (1e-8) diagonal floor of the band retry matrix, as before.

## Pending refactorings

- `solve()` is a 1000-line function mixing ordering, workspace setup, fill
  and five solver routes; the fill loops should become one helper per
  storage format.
