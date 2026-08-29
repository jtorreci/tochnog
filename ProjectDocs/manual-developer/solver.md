# solver (the global solver family) — developer notes

Sprint 12 lot 3. The Professional's global solver family (manual
6.1047-6.1056) as PLAIN (non-indexed) records - unlike their
`control_solver_*` counterparts of Sprint 9, which are CONTROL-class
(indexed) records.

## Implementation

- **Enums**: `SOLVER`, `SOLVER_BICG_ERROR`, `SOLVER_BICG_RESTART`,
  `SOLVER_BICG_STOP`, `SOLVER_MATRIX_SAVE`, `SOLVER_MATRIX_SYMMETRIC`,
  `SOLVER_PARDISO_ORDERING`, `SOLVER_PARDISO_OUT_OF_CORE`,
  `SOLVER_PARDISO_PROCESSORS`, `SOLVER_PARDISO_PROCESSORS_MAXIMUM` in
  `tochnog.h`/`tochnog-mod.h` (synced, same order). Registered in
  `database.cc` with `no_index = 1` (plain records; the Professional
  synopsis has NO record index).

- **`solver` (global type override)**: read LAST in the three
  consumption points of the solver type, so it overwrites both the
  default and the per-control record:
  - `top.cc` (main solve dispatch, after `OPTIONS_SOLVER` and
    `CONTROL_OPTIONS_SOLVER` - note the GNU's own `options_solver`
    global record has the OPPOSITE precedence: the control record
    overwrites IT; the Professional's `solver` is the other way
    around, hence the extra read),
  - `elem.cc` (element matrix assembly branch),
  - `dof.cc` (dof update branch).
  Proof: `tslv_override` - with `control_solver 0
  -matrix_iterative_bicg` and `solver -matrix_superlu` the SuperLU
  banner appears (the direct solver ran).

- **`solver_bicg_error`**: read in `so_bicg.cc` right after
  `CONTROL_OPTIONS_SOLVER_BICG_ERROR`, overwriting it.

- **`solver_bicg_stop`**: read in `so_bicg.cc` after the
  `CONTROL_SOLVER_BICG_STOP` block, overwriting it.

- **`solver_matrix_symmetric -yes`**: in `so_bicg.cc`, replaces the
  per-solve measurement `solve_iterative_bicg_symmetric()` with the
  user assertion (plain CG). The measurement still runs (once) to
  print an unconditional WARNING when it disagrees - the Professional
  symmetrizes the matrix in that case, the GNU runs CG on the
  as-assembled matrix, so the warning tells the user CG may not
  converge (test `tslv_bsym`: the beam, measured non-symmetric by the
  dtime factor in the translation-rotation blocks, prints the warning
  and survives with `solver_bicg_stop -no`). `-no`/absent: measured
  behavior unchanged (CG when symmetric, honest Bi-CG when not).

- **PARTIALS** (registered, accepted, behavior-neutral - test
  `tslv_part` verifies identical results with all of them present):
  `solver_bicg_restart` (the GNU bicg has no restart),
  `solver_matrix_save` (decomposition caching is direct-solver
  functionality; `-always` is pardiso-only even in the Professional),
  the four `solver_pardiso_*` (PARDISO is not compiled in the GNU;
  `solver_pardiso_ordering` is registered INTEGER - the Professional
  value is an ordering name).

## Gotchas

- The records are PLAIN: `no_index = 1`. Giving them a control index
  (`solver 0 -matrix_superlu`) misparses (0 would be read as the
  type).
- `zip`-style system() calls: see `zip.md` developer notes for the
  nullglob trap.
