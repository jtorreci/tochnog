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

## Measurements (omen, one input at a time, 12 GB cap, timeout 600 s)

Professional `other` inputs; old = `225d536`, new = `ab99338` (A+B+C);
Professional 02-08-2026 from the F4 T7 investigation (same script).
Peak RSS from `/usr/bin/time -v` (new `force11` / `foundation...`: VmHWM
sampled over 140 s, the runs hit the timeout).

| input | old rc / RSS / wall | new rc / RSS / wall | Pro rc / RSS / wall |
|---|---|---|---|
| `force11` | 134 bad_alloc / 7.8 GB / 43 s | 124 timeout / 0.27 GB / >600 s | 0 / 0.27 GB / 13 s |
| `large1` | 134 bad_alloc / 0.65 GB / 25 s | 0 / 0.65 GB / 58 s | 0 / 1.3 GB / 9 s |
| `ground19_water_under_dam` | 0 / 5.2 GB / 20 s | 0 / 0.12 GB / 18 s | 0 / 0.17 GB / 3 s |
| `mesh_activate1` | 134 bad_alloc / 11.9 GB / 279 s | 1 (target, node numbering) / 0.035 GB / 29 s | 0 / 0.07 GB / 26 s |
| `foundation_consolidation_settlement` | 134 bad_alloc / 0.19 GB / 7 s | 124 timeout / 1.5 GB / >600 s | 124 timeout / - / >600 s |

Results vs the Professional: `ground19` meets its target (rc 0) as the
Professional; `mesh_activate1` 202 nodes and integrated `veliy` equal to
5e-10 (its target names a node number of the Professional macro
numbering; final `vely` differs, BACKLOG C46); `large1` has no target
and no database. Open: `force11` (one 3D solve) and
`foundation_consolidation_settlement` (Bi-CG breaks down, band LU retry
on every solve) are now limited by the solver speed, not memory
(BACKLOG C43).

Corpus wall time (omen, 6 jobs, tests with rc 0 on both binaries), old
-> new: Pro other 590.2 -> 552.8 s (250 tests), sfnet 182.9 -> 176.2 s
(173), test-2026 35.6 -> 24.5 s (99, `mesh_macro_timestep1` 11 -> 0.9 s),
validation-suite 29.9 -> 31.3 s (233, noise); total 838.6 -> 784.8 s
(-6 %). Largest gains: `force5` 8.0 -> 3.0 s, `force9` 51 -> 46 s,
`taylor3` 92 -> 84 s, `large2` 55 -> 50 s.

## External features / libraries

- LAPACK `dgbsv_` (band LU).

## Hardcoded parameters

- `EPS_MAT` (1e-8) diagonal floor of the band retry matrix, as before.

## Pending refactorings

- `solve()` is a 1000-line function mixing ordering, workspace setup, fill
  and five solver routes; the fill loops should become one helper per
  storage format.

## F4 T9 (2026-10-10): solver speed

Changes (BACKLOG C43, owner decisions 2026-10-10): `-matrix_pardiso` is
solved with SuperLU with a Bi-CG fallback ([matrix_pardiso](matrix_pardiso.md),
`3349515`); the CG/Bi-CG products run on an assembled CSR matrix
([solver_bicg_csr_operator](solver_bicg_csr_operator.md), `b7009dd`); the
nodal stress recovery factors are tabulated (`acfac40`); two SuperLU /
interface memory defects fixed (`681dc83`, `789f8bf`). The GNU default
stays the Bi-CG.

Measurements (omen, one input at a time, 12 GB cap, timeout 900 s,
`/usr/bin/time -v`; old = `c62d28f`, new = `b7009dd`, md5 e84e554c;
Professional 02-08-2026 with the same script). A run killed by the
timeout has no `/usr/bin/time` record: its peak RSS is the VmHWM sampled
in F4 T7 (old) or locally (new `foundation...`, 4 GB cap, 900 s).

| input | old rc / RSS / wall | new rc / RSS / wall | Pro rc / RSS / wall |
|---|---|---|---|
| `force11` (`solver -matrix_pardiso`) | 124 timeout / 0.27 GB (F4 T7) / > 900 s (Bi-CG; > 1 h estimated) | 1 (target C47) / 0.41 GB / 14.9 s (SuperLU) | 0 / 0.27 GB / 11.9 s (PARDISO) |
| `large1` (default Bi-CG) | 0 / 0.65 GB / 60.0 s | 0 / 0.72 GB / 23.9 s | 0 / 1.31 GB / 9.6 s |
| `foundation_consolidation_settlement` | 124 timeout / 1.5 GB (F4 T7) / > 900 s | 124 timeout / 1.56 GB / > 900 s | 124 timeout / - / > 600 s (F4 T7) |

`large1`: -60 % wall (recovery table 37 -> 27 s locally, the rest from
the CSR products) for +70 MB (the CSR matrix). `force11` with the CSR
Bi-CG (`-matrix_iterative_bicg`, local): 311040 iterations (the maximum)
in ~640 s without converging (the Professional Bi-CG does not converge
either), then the band LU retry: same `mom1y_sig` 5.80223 as SuperLU.
`foundation_consolidation_settlement` stays Bi-CG bound (breakdown and
band LU retry on every solve); the Professional times out as well.

Corpus wall time (omen, 6 jobs, 4 GB cap, timeout 120, sum of the test
times of tests finished on both), old -> new: Pro other 794 -> 525 s
(317 tests), sfnet 317 -> 220 s (195), test-2026 39 -> 33 s (113),
validation-suite 31 -> 24 s (242); total 1180 -> 803 s (-32 %). Per
step: recovery table (old -> b1) Pro other 794 -> 685 s, sfnet 317 ->
289 s; CSR (b3 -> b4) Pro other 696 -> 545 s, sfnet 288 -> 220 s. OK
counts: Pro other 273 -> 276 (`force11`, `ground11_nonsaturated`,
`ground11_phreatic_level` finish), test-2026 110 -> 111
(`superlu_solve1`), sfnet 192 and validation-suite 232 unchanged.
