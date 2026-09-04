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

## 2026-09-04: direct SuperLU route fixed (was segfaulting) — `c2cbaf7`

`so.cc` + `so_suplu.c`. Two stacked root causes for the
`solver -matrix_superlu` segfault on penalty interface models
(patch1, interface_bar2_hex8):

1. **Dead sparse assembly since 5b7bb8b.** That commit started to
   allocate the band workspace `mat` in EVERY route (needed by the
   direct-LU retry after a Bi-CG failure). The assembly branch was
   `if ( mat ) { band fill } else { assert(petsc||superlu); sparse fill }`,
   so with `mat` always allocated the sparse fill NEVER ran and
   `solve_superlu()` received an EMPTY matrix (`superlu_nnz == 0`,
   verified with a temporary dump) -> crash inside SuperLU. The sparse
   branch now runs whenever the active solver is SuperLU/PETSc
   (`if ( mat && !petsc_solver && !superlu_solver )`).

2. **Wrong-stride membership search in the sparse fill.**
   `array_member( (long int *) inz[a], b, nnz[a], index )` scanned the
   `int` row buffers as `long int` (8-byte stride). For columns with
   more than two entries the 64-bit composite read of two consecutive
   `int`s never matched the searched row index, so repeated (row,col)
   pairs were NOT merged: duplicate CSC entries, values accumulated
   into wrong slots and missing diagonals resulted (measured on
   bar2_hex8: nnz dropped 437 -> 376 once merged, 47/60 columns had
   duplicates, column 1 lost its diagonal). SuperLU's colamd then
   crashed on the malformed structure. Replaced by the int-correct
   helper `sparse_list_add()` (append-or-accumulate with growth).

Additional changes:

- `mat`/`ipiv` (band) are freed in every route at the end of
  `solve()` (previously only with `band_solver`: per-call leak on the
  Bi-CG and SuperLU routes since 5b7bb8b).
- `so_suplu.c`: `solve_superlu()` now returns `0` when `dgssv` fails
  (`info != 0`, e.g. a singular matrix) so `so.cc` aborts with an
  honest error instead of continuing with a garbage solution;
  `dQuerySpace()` runs before the `mem_usage` print (it was read
  uninitialized); `StatFree()` added after `dgssv`.

### Verification

- patch1 with `solver -matrix_superlu`: rc=0, four real factorizations
  ("Solution Found" banners), interface sigma_n identical to the band
  LU route (element 5: 1431.07742001 vs 1431.07742002).
- interface_bar2_hex8 with `solver -matrix_superlu`: no crash, honest
  `dgssv` INFO=14 (the system is genuinely singular - 18 zero modes -
  a mesh/extrusion bug of that test, unrelated to the solver).
- Reduced suite 16/16 (`tslv_override` exercises the SuperLU banner).
- Blast radius (114 solver/interface tests): 59 PASS, unchanged vs
  HEAD; force16/force17/mpc8/mpc9 segfault at HEAD as well
  (pre-existing, MPC + stress-dof family, unrelated to this change).
