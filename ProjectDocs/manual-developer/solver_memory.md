# Solver memory — developer notes (F4 T7)

## Where it is implemented

- `so.cc`, function `solve()`: workspace allocation (after the eigen
  checks), the band/sparse fill loops and the Bi-CG failure branch.

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

## External features / libraries

- LAPACK `dgbsv_` (band LU).

## Hardcoded parameters

None.

## Pending refactorings

- `solve()` is a 1000-line function mixing ordering, workspace setup, fill
  and five solver routes; the fill loops should become one helper per
  storage format.
