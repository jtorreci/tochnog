# control_print_vtk_empty

## Implementación

- Keyword `control_print_vtk_empty` (data_class CONTROL, type INTEGER,
  data_length 1) registered in `database.cc` with
  `data_required = CONTROL_PRINT_VTK`.
- New enum `CONTROL_PRINT_VTK_EMPTY` in `tochnog.h` / `tochnog-mod.h`
  (kept in sync).
- Helper `vtk_element_is_empty(element)` in `print_vt.cc`: reads
  `ELEMENT_EMPTY` (VERSION_PRINT, `GET_IF_EXISTS`) and returns true when
  the value is `-YES` — the pattern of `print_g5.cc` (an element is
  present when the record is `-NO` or `-FRONT`).
- With `-no` the empty elements are skipped in the three element passes:
  the `CELLS` length/count pass (also computes `ncell`), the `CELLS`
  write and the `CELL_TYPES` write. The header counts (`CELLS ncell`,
  `CELL_TYPES ncell`) stay consistent with the written cells.

## Qué es un elemento vacío

- `ELEMENT_EMPTY` is COMPUTED by the solver (elem.cc), not user input:
  with `materi_diffusion` an element is `-YES` when NONE of its nodes
  has diffusion >= `EPS_MATERI_DIFFUSION_MINIMUM` (0.5); with
  `materi_density` when NONE has density >=
  `EPS_MATERI_DENSITY_MINIMUM` (1e-9). Empty elements are skipped by
  the element loop (they have no stiffness/results).
- The record is written to VERSION_NEW in the element loop and reaches
  VERSION_NORMAL through the NEW→NORMAL copy before `step_close`; it
  has `version_all=1` so it is available in VERSION_PRINT.

## Decisiones

- The "illegal element type" abort for unsupported element names is
  UNCHANGED: with `-no` the empty elements are skipped BEFORE the type
  check, so an empty element with an exotic type no longer aborts (it
  is simply not written); with the default `-yes` the behavior is
  byte-identical to before. No legacy test regressed (155/155).
- `CELLS` connectivity uses node indices, so skipping an element does
  NOT renumber the points (only the cell count changes).
- Verification: `vtk_empty1` test — element 2 empty via node density 0;
  default `tn35.vtk` has `CELLS 2`, `-no` `tn36.vtk` has `CELLS 1 5`
  and `CELL_TYPES 1`. The model still solves (nodes without elements
  are handled by the solver, pattern of the deletion tests).

## Pendiente

- Nothing. The empty-element detection is the solver's own
  `ELEMENT_EMPTY` record; there is no user-facing keyword to force an
  element empty.
