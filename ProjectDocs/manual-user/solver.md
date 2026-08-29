# solver (the global solver family)

## Description

The `solver` data record sets the solver type for the ENTIRE
calculation. Unlike `control_solver`, which only holds for the time
steps of its own control index, a specified `solver` OVERWRITES every
`control_solver` record (manual Professional 6.1047).

The record is plain (no index):

```
solver -matrix_superlu
```

Valid types are the same as for `control_solver`:
`-matrix_iterative_bicg` (default), `-matrix_superlu`,
`-matrix_lapack` (band), and the PETSC variants when compiled in.

When using the bicg solver, consider also setting
`solver_matrix_symmetric -yes` to speed up the solve.

## The global solver_* family (manual Professional 6.1048-6.1056)

These plain records are the global counterparts of the per-control
`control_solver_*` records of Sprint 9. When both exist, the global
record overwrites the per-control one.

| record | meaning | GNU status |
|---|---|---|
| `solver_bicg_error <error>` | termination error ratio of the bicg iterations (6.1048) | wired (overwrites `control_solver_bicg_error`) |
| `solver_bicg_restart <n>` | number of bicg restarts (6.1049) | PARTIAL: registered and accepted; the GNU bicg has no restart mechanism |
| `solver_bicg_stop <sw>` | stop the calculation when bicg does not converge (6.1050) | wired (overwrites `control_solver_bicg_stop`; `-no` continues with a warning) |
| `solver_matrix_save <sw>` | save and reuse the decomposed matrix (6.1051) | PARTIAL: registered and accepted; decomposition caching is direct-solver functionality not present in the GNU. `-always` is pardiso-only in the Professional |
| `solver_matrix_symmetric <sw>` | treat the matrix as symmetric (6.1052) | wired, see below |
| `solver_pardiso_ordering <o>` | pardiso ordering (6.1053) | PARTIAL: registered; PARDISO is not compiled in the GNU |
| `solver_pardiso_out_of_core <sw>` | pardiso out-of-core (6.1054) | PARTIAL, idem |
| `solver_pardiso_processors <n>` | pardiso processor count (6.1055) | PARTIAL, idem |
| `solver_pardiso_processors_maximum <sw>` | pardiso max processors (6.1056) | PARTIAL, idem |

## solver_matrix_symmetric

With `-yes` the user asserts that the system is symmetric: the
per-solve symmetry measurement is bypassed and the symmetric solver
(plain CG) is used. The GNU measures the symmetry of every solve
(measured symmetric -> CG, non-symmetric -> honest Bi-CG); this record
removes the measurement from the equation. The Professional
symmetrizes the matrix when needed; the GNU runs CG on the
as-assembled matrix, so a WARNING is printed when the measurement
disagrees with the user's assertion:

```
Warning: solver_matrix_symmetric -yes but the measured system is NOT
symmetric - CG may not converge on it
```

Combine with `solver_bicg_stop -no` to keep the calculation alive in
that situation.

## Example

```
control_solver 0 -matrix_iterative_bicg
solver -matrix_superlu          # wins over the control record
solver_matrix_symmetric -yes    # (only sensible with an iterative type)
solver_bicg_stop -no
```

The `tslv_override` test of the validation suite proves the override:
with `control_solver` asking for bicg and `solver -matrix_superlu`,
the SuperLU banner appears in the output (the direct solver ran).
