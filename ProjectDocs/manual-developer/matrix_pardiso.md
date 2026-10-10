# matrix_pardiso — developer notes (F4 T9, 2026-10-10)

## Where

`so.cc` `solve()`:

- the `task==-MATRIX_PARDISO` branch sets `superlu_solver` and
  `pardiso_superlu` (Bi-CG when SuperLU is not compiled in, or after an
  earlier failure: `solve_pardiso_superlu_failed`);
- the solver call section now runs SuperLU first (for `-matrix_superlu`
  and `-matrix_pardiso`); on a failed `-matrix_pardiso` factorization it
  restores the right-hand side from `solve_b_temp`, zeroes `solve_x`,
  sets `bicg_solver` and `solve_pardiso_superlu_failed`, and the
  following Bi-CG block (with the band LU retry `solve_band_retry`) solves
  the same system;
- `solve_pardiso_memory_guard()` / `solve_proc_kb()`: the memory guard.

Owner decision 2026-10-10 (BACKLOG C43): `-matrix_pardiso` -> SuperLU
with Bi-CG fallback; the GNU default stays the Bi-CG (SuperLU does not
scale in 3D: the Professional `large1` 30^3 needs > 4 GB); C12/C13
(`bicg_stop`) semantics unchanged (they act inside
`solve_iterative_bicg()`).

## Memory guard

`dgssv` keeps its `GlobalLU_t` local, so when the factorization runs out
of memory (`info > n`) the L/U workspace it allocated is never released.
Before the factorization `solve_pardiso_memory_guard()` lowers the soft
`RLIMIT_AS` to `VmSize + PARDISO_SUPERLU_MEMORY_FRACTION * room`, with
`room = min(soft limit - VmSize, MemAvailable)` (from `/proc/self/status`
and `/proc/meminfo`); the limit is restored right after `dgssv`. SuperLU
then fails cleanly at that budget and the remaining 30 % of the room is
left for the Bi-CG. A failed factorization disables SuperLU for the rest
of the run (each failure leaks its budget). Without `/proc` or without
room no guard is set. `TN_PARDISO_SUPERLU_MB` replaces the budget
(diagnostics, test `pardiso_bicg_fallback1`).

Residual risk: SuperLU aborts the program (`ABORT`) when some of its
small allocations fail (`get_perm_c`/`colamd` work arrays, the `dgstrs`
work vector); with the 70 % budget these are far smaller than the matrix
the GNU already holds.

## Hardcoded parameters

- `PARDISO_SUPERLU_MEMORY_FRACTION` 0.7.
- SuperLU options: `set_default_options` (COLAMD column ordering,
  partial pivoting threshold 1.0), as `-matrix_superlu`.

## Verification

- `pardiso_superlu1` (13x7x7 hex8, analytic): SuperLU solves it ("Solution
  Found"); old GNU: Bi-CG ("falling back"). Results equal to the Bi-CG to
  the Bi-CG tolerance (disy -0.00120000000011 vs -0.001200000002).
- `pardiso_bicg_fallback1` (25x13x13 hex8) with `TN_PARDISO_SUPERLU_MB=4`:
  `dgstrf` info > n, fallback to the Bi-CG, targets met.
- Professional `large1` with `solver -matrix_pardiso` under a 4 GB cap:
  SuperLU fails cleanly (info > n), the Bi-CG finishes, rc 0, 47 s,
  2.54 GB peak (no crash; before the so_suplu.c fix a failed
  factorization crashed).
- Professional `force11` (`solver -matrix_pardiso`): 10 s, 0.41 GB
  (old: Bi-CG > 1 h). The run ends with rc 1 on its `mom1y_sig` target
  (5.80223 vs 5.76474; Professional 5.76460): a 0.65 % model difference
  of `post_calcul -materi_stress -force` in 3D exposed now that the solve
  finishes (BACKLOG C47). Cause found in F4 T12 (2026-10-11): the
  quad9 -> hex27 extrusion disconnected the elements at the inner
  planes (`control_mesh_extrude`); after the fix rc 0, 5.76461.
- Corpus A/B and timings: SEGUIMIENTO F4 T9 row and
  `solver_memory.md`.

## Pending refactorings

- The SuperLU sparse fill (`sparse_list_add`, linear scan per column) and
  the separate CSC copy could reuse the CSR pattern of `so_bicg.cc`.
- The repository compiles `so_suplu.c` against SuperLU 4.1 headers while
  linking SuperLU 6.0.1 (see `solver.md`); the SuperLU 6.0.1 headers
  would allow `dgstrf` with a user work array (`lwork > 0`) instead of
  the address-space guard.
