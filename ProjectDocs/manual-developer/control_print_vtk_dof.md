# control_print_vtk_dof

## Implementación

- **Filter**: applied in `print_vtk()` in `print_vt.cc`. When
  `control_print_vtk_dof` is active for the control record, only the
  solution fields whose `dof_type` matches one of the listed initialisation
  names are written. `-none` writes no field.
- **Keyword**: `control_print_vtk_dof` (data_class CONTROL, type INTEGER,
  data_length DATA_ITEM_SIZE, fixed_length 0) registered in `database.cc`,
  with `data_required = CONTROL_PRINT_VTK`.
- **New enum**: `CONTROL_PRINT_VTK_DOF` in `tochnog.h` / `tochnog-mod.h`
  (kept in sync), right after `CONTROL_PRINT_VTK`.
- **Three write loops respect the filter**:
  1. Primary unknowns (SCALARS/VECTORS/TENSORS headers + node data).
  2. Vector/tensor component scalars.
  3. Derived magnitudes (von Mises, Tresca, principal stresses) — only
     written when the filter (if any) includes `-materi_stress`.

## Diseño / decisiones

- `print_field` is recomputed per `iuknwn` from `dof_type[iuknwn]`
  compared against the listed names. `nval` is always computed (even when
  filtered out) so `ipuknwn += nval` advances correctly.
- The derived-magnitude block checks `-MATERI_STRESS` explicitly; without
  a filter it is always written (default behavior unchanged).

## Detalles

- The filter matches on `dof_type` (e.g. `-MATERI_STRESS`), not on the
  scalar/vector component labels.
- `-none` is detected via `vtk_dof[0]==-NONE`.

## Pendiente

- `control_print_vtk_coord`, `_empty`, `_node_method`, `_other` variants
  are not implemented (the Professional family has more sub-options).
