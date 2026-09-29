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

`-matrix_superlu` runs the direct sparse SuperLU LU factorization
(2026-09-04: the route was segfaulting; fixed in `c2cbaf7`). It is the
direct-solver alternative to `-matrix_pardiso` (PARDISO itself is not
compiled in the GNU - that type falls back to the iterative Bi-CG with
an automatic retry with the direct band LU when the Bi-CG fails).

## The global solver_* family (manual Professional 6.1048-6.1056)

These plain records are the global counterparts of the per-control
`control_solver_*` records of Sprint 9. When both exist, the global
record overwrites the per-control one, except `solver_bicg_stop`, which
the control record overrules (manual 6.376, see the table).

| record | meaning | GNU status |
|---|---|---|
| `solver_bicg_error <error>` | termination error ratio of the bicg iterations (6.1048) | wired (overwrites `control_solver_bicg_error`) |
| `solver_bicg_restart <n>` | number of bicg restarts (6.1049) | PARTIAL: registered and accepted; the GNU bicg has no restart mechanism |
| `solver_bicg_stop <sw>` | stop the calculation when bicg does not converge (6.1050) | wired (`control_solver_bicg_stop` overrules it: the index-0 record, then the record of the timestep's own index; `-no` continues with a warning) |
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

## Differences with the Professional version

- `control_solver_bicg_stop` / `solver_bicg_stop` precedence. The GNU
  reads, last one wins: `solver_bicg_stop`, then
  `control_solver_bicg_stop 0 ...`, then `control_solver_bicg_stop` with
  the control index of the running timestep. The Professional
  (02-08-2026, measured) reads only the index-0 control record, for every
  timestep, and ignores a record with the timestep's own non-zero index,
  although the manual (6.376) documents a per-timestep index. The GNU
  follows the manual here (project decision 2026-09-29): inputs with
  only `control_solver_bicg_stop 0 ...` and/or `solver_bicg_stop` behave
  as in the Professional; an input with a record at a non-zero timestep
  index behaves differently. Example, a timestep with index 20:
  `control_solver_bicg_stop 20 -no` continues in the GNU and stops in
  the Professional; `control_solver_bicg_stop 0 -no` plus
  `control_solver_bicg_stop 20 -yes` stops in the GNU and continues in
  the Professional.
