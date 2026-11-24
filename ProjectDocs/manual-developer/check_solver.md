# check_solver

## Archivos y funciones

- `so.cc` → `solve( long int task )` (`so.cc:58`) — the check lives inside the
  `if (band_solver)` block (`so.cc:639`) just before the `dgbsv_` call
  (`so.cc:656`). It loops over the `solve_nlocal` equations and compares
  `scalar_dabs( mat[kl+ku + idiag*ldmat] )` against `check_solver_eps`
  (`so.cc:644-655`).
- `initia.cc:96` — global `double check_solver_eps=0.;` (default: check off).
- `top.cc:130` — reads the keyword once at startup with
  `db( CHECK_SOLVER, 0, idum, &check_solver_eps, ldum, VERSION_NORMAL, GET_IF_EXISTS )`.
  Note the value arrives through the `dval` argument because the type is
  `DOUBLE_PRECISION`, not `ival`.
- `database.cc:297-300` — keyword registration:
  `strcpy(name[CHECK_SOLVER],"check_solver")`, `type = DOUBLE_PRECISION`,
  `data_length = 1`, `no_index[CHECK_SOLVER] = 1`.
- Enum `CHECK_SOLVER` in `tochnog.h:169` / `tochnog-mod.h:162` (must stay in sync).

## Detalles de implementación

- The check only runs when `band_solver` is set, i.e. with
  `options_solver -matrix_lapack` (`so.cc:163-167`).
- The diagonal of equation `i` in LAPACK band storage is at index
  `kl+ku + i*ldmat` (row offset `kl`, column offset `ku` on band row `i`),
  matching the layout expected by `dgbsv_`.
- Guarded by `check_solver_eps>0.`: a non-positive value disables the check.
- On a hit it prints the equation number (`idiag+1`), the diagonal value and
  `eps` to `std::cout`, plus the hint
  `This normally indicates a problem in the input file.` It does not abort the
  run; the solver proceeds with the `dgbsv_` call regardless.

## Dependencias externas

None beyond LAPACK itself: the block sits inside the existing `dgbsv_` call
path, so the feature only compiles/works when LAPACK is linked. Uses
`scalar_dabs` and the `band_solver` local inside `solve()`.

## Parámetros hardcodeados / refactorizaciones pendientes

- The warning prints `mat[kl + idiag*ldmat]` while the check compares
  `mat[kl+ku + idiag*ldmat]`; on `ku>0` the printed value is NOT the tested
  diagonal element — the printed index should be `kl+ku` too.
- Output goes to `std::cout`; routing it through `pri()` / `tn.log` would be
  consistent with the rest of the solver error reporting.
- `check_element_shape` and `check_memory` remain unimplemented in the same
  `check_*` family.
