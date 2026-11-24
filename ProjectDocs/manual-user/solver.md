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
`-matrix_pardiso` (solved with SuperLU, see below), `-matrix_lapack`
(band), and the PETSC variants when compiled in.

When using the bicg solver, consider also setting
`solver_matrix_symmetric -yes` to speed up the solve.

`-matrix_superlu` runs the direct sparse SuperLU LU factorization
(2026-09-04: the route was segfaulting; fixed in `c2cbaf7`; 2026-10-10,
F4 T9: a second crash after the first solve of larger models - a read of
the factors after they were released - is fixed, test `superlu_solve1`). PARDISO itself
is not compiled in the GNU: since F4 T9 (2026-10-10) `-matrix_pardiso` is
solved with SuperLU, with a fallback to the iterative Bi-CG (and its
direct band LU retry) when the factorization fails - see
[matrix_pardiso](matrix_pardiso.md). The GNU default stays the Bi-CG.

### Stopping criterion of the iterative solver (F4 T6f C77, 2026-10-15)

The iterative solver (CG for a symmetric system, Bi-CG otherwise) stops
when the residual of the current iterate, relative to the right-hand
side, is below the tolerance:

```
|b - A x|^2 / |b|^2  <  bicg_error      (Jacobi-scaled norm, default 1e-14)
```

`bicg_error` is `control_solver_bicg_error` / `solver_bicg_error` (the
termination error ratio, manual 6.374 / 6.1048), so the default is a
relative residual of 1e-7. The stop does not depend on the scale of the
load: a tiny correction solve is solved as accurately as the first one.

`control_options_solver_bicg_error_minimum <value>` (and the global
`options_solver_bicg_error_minimum`, GNU records) add an ABSOLUTE floor:
the solver also stops when `|b - A x|^2 < value` (Jacobi-scaled). Without
the record there is no absolute floor. Before 2026-10-15 a floor of 1e-16
was always active; a solve with a small right-hand side (a Newton
correction, a flow near its steady state) then stopped with a relative
residual up to 0.13.

Two more stops protect the iteration:

- the residual cannot drop below the round-off level of the computer,
  about `1e-14 |A| |x|`; the iteration stops there;
- the iteration breaks down when a scalar product of the method becomes
  zero RELATIVE to the size of its vectors.

The iterate of such a stop is kept when its relative residual
`|b - A x|/|b|` is at most 1e-6. Otherwise the system is re-solved with a
direct solver: SuperLU first, then the direct band LU when SuperLU fails
and the band matrix fits in memory. The log shows one line per such
solve:

```
Warning: iterative solver breakdown, system solved with the direct SuperLU
```

(or `... with the direct band LU`). This is information, not an error in
the input. When the iterative solver does not converge within its
iteration limit, the log shows once "Warning: the iterative solver did
not converge - retrying with a direct solver". When no direct solver can
solve the system either, the run stops with the full "Error: the
iterative solver broke down ... Probably you have a wrong input file"
block, or with "Error: the direct band LU retry does not fit in memory"
when the band matrix is too large (before 2026-10-15 a large model
crashed here with std::bad_alloc). With `solver_bicg_stop -no` the
partial iterative solution is kept, with a warning.

Example (validation_1 channel flow, `convection_apply -yes`): before
2026-10-15 the GNU printed 599 band-LU warnings (false breakdowns of the
old absolute test, F4 T6f C74); now every solve converges. With the
inflow scaled by 1e-5 the old floor gave velx 9.77556e-7, now 9.775171e-7
= SuperLU (test `bicg_relative_stop1`).

## The global solver_* family (manual Professional 6.1048-6.1056)

These plain records are the global counterparts of the per-control
`control_solver_*` records of Sprint 9. When both exist, the global
record overwrites the per-control one, except `solver_bicg_stop`, which
the control record overrules (manual 6.376, see the table).

| record | meaning | GNU status |
|---|---|---|
| `solver_bicg_error <error>` | termination error ratio of the bicg iterations (6.1048) | wired (overwrites `control_solver_bicg_error`): relative stop `|r|^2/|b|^2 < error`, see above |
| `solver_bicg_restart <n>` | number of bicg restarts (6.1049) | PARTIAL: registered and accepted; the GNU bicg has no restart mechanism |
| `solver_bicg_stop <sw>` | stop the calculation when bicg does not converge (6.1050) | wired (`control_solver_bicg_stop` overrules it: the index-0 record, then the record of the timestep's own index; `-no` continues with a warning) |
| `solver_matrix_save <sw>` | save and reuse the decomposed matrix (6.1051) | PARTIAL: registered and accepted; decomposition caching is direct-solver functionality not present in the GNU. `-always` is pardiso-only in the Professional |
| `solver_matrix_symmetric <sw>` | treat the matrix as symmetric (6.1052) | wired, see below |
| `solver_pardiso_ordering <o>` | pardiso ordering (6.1053) | PARTIAL: registered; PARDISO is not compiled in the GNU |
| `solver_pardiso_out_of_core <sw>` | pardiso out-of-core (6.1054) | PARTIAL, idem |
| `solver_pardiso_processors <n>` | pardiso processor count (6.1055) | PARTIAL, idem |
| `solver_pardiso_processors_maximum <sw>` | pardiso max processors (6.1056) | PARTIAL, idem |

## solver_matrix_symmetric

With `-yes` the symmetric solver (plain CG) is used when the system
measures symmetric. The GNU measures the symmetry of every solve
(measured symmetric -> CG, non-symmetric -> honest Bi-CG) and does NOT
symmetrize the matrix: when the measurement disagrees with `-yes`, it
runs Bi-CG on the as-assembled matrix and prints a WARNING:

```
Warning: solver_matrix_symmetric -yes but the measured system is NOT
symmetric - running Bi-CG on the as-assembled matrix instead of CG (the
Professional symmetrizes the matrix in this case)
```

The GNU result therefore does not depend on the switch. The
Professional symmetrizes the matrix; on coupled consolidation (u-p)
systems this changes its answer (deliberate divergence K24, see
[solver_matrix_symmetric](solver_matrix_symmetric.md)).

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
- `solver_matrix_symmetric -yes` on a non-symmetric system (coupled
  consolidation): the Professional symmetrizes the matrix, the GNU
  solves the as-assembled one (K24, deliberate divergence, project
  decision 2026-10-11; details in
  [solver_matrix_symmetric](solver_matrix_symmetric.md)).

- **Iterative solver stop** (C74, C77; DIVERGENCES C77, class A). The
  Professional default solver is the direct PARDISO solver. Its iterative
  `solver -matrix_iterative_bicg` stops on an absolute convergence floor:
  `validation_1` with the inflow scaled by 1e-5 gives velx 8.9888e-7
  against PARDISO 9.77518e-7. The GNU iterative solver stops on the
  relative residual (owner decision 2026-10-15): 9.775171e-7, equal to
  SuperLU. The Intel MKL PARDISO solver is not integrated in the GNU
  (decision 2026-10-15, SEGUIMIENTO Descarte); `-matrix_pardiso` is
  solved with SuperLU. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **Bi-CG breakdown on a singular system** (BACKLOG C12, closed
  2026-10-17, no divergence): re-measured on an unsupported elasti7 mesh,
  the Professional 02-08-2026 finishes with rc 0 and a `.dbs` (no
  message); the GNU prints `Warning: iterative solver breakdown, system
  solved with the direct SuperLU` (C77 fallback) and also finishes with
  rc 0. Neither code detects the singular model.
