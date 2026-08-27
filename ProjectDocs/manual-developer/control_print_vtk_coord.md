# control_print_vtk_coord

## Implementación

- Keyword `control_print_vtk_coord` (data_class CONTROL, type INTEGER,
  data_length 1) registered in `database.cc` with
  `data_required = CONTROL_PRINT_VTK`.
- New enum `CONTROL_PRINT_VTK_COORD` in `tochnog.h` / `tochnog-mod.h`
  (kept in sync), right after `CONTROL_PRINT_VTK`.
- Read in `print_vtk()` (`print_vt.cc`) with `GET_IF_EXISTS`; default
  `-YES` when the record is absent. The `POINTS` block (the only place
  where coordinates appear) is written only when the switch is not `-NO`.

## Decisión (documentada, manual 6.340)

- **Evidence from the code**: before this feature, `print_vt.cc` wrote
  the coordinates ONLY in the `POINTS` block (the mesh section). There
  are NO coordinate data fields in `POINT_DATA` — unlike gmsh (manual
  6.316 describes `node_0`/`node_1`/`node_2` scalars, which the GNU
  `print_gm.cc` does not implement either).
- **Interpretation**: "the coordinates of nodes is not plotted in vtk"
  = the `POINTS` block is omitted. The alternative (adding coordinate
  data fields like gmsh) was rejected because it would change the
  default output of an existing feature and the GNU gmsh writer does not
  write coordinate fields either.
- **Known limitation**: a VTK UNSTRUCTURED_GRID without `POINTS` is not
  a valid dataset for a post-processor (cells reference non-existent
  point indices). The switch is a file-size/debug tool; documented in
  the user manual.
- The "illegal element type" aborts and every other block are untouched:
  with `-no` only the coordinate section disappears.

## Detalles

- The `POINTS` count is still `max_node+1` (the node numbering is
  untouched by the switch).
- Verification: `vtk_coord1` test — A/B `tn30.vtk` (default, has
  `POINTS`) vs `tn31.vtk` (`-no`, no `POINTS` but still `CELLS`).

## Pendiente

- Nothing (switch fully implemented). The gmsh-style coordinate data
  fields (`node_0`/`node_1`/`node_2`) remain unimplemented (out of
  scope, would change default output).
