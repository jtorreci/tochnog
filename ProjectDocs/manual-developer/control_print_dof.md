# control_print_dof

## Implementación

- **Output**: `print_dof()` in `print_hi.cc` (same file as
  `print_history` / `print_history_smooth`). Invoked from the control loop
  in `top.cc`:
  ```
  if ( db_active_index( CONTROL_PRINT_DOF, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_DOF, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
    if ( ival[0]!=-NO ) print_dof( icontrol, ival[0] );
  }
  ```
- **Keyword**: `control_print_dof` (data_class CONTROL, data_length 1,
  type INTEGER) registered in `database.cc`. `ival[0]` holds the switch.
- **New enum**: `CONTROL_PRINT_DOF` in `tochnog.h` / `tochnog-mod.h`
  (kept in sync).
- **Layout**: for each time step (each call appends), for each dof
  component (in `node_dof` order via `dof_scal_vec_mat`), for each node,
  write `x y z <value>`. Matrices (stress) use `stress_indx(kk,ll)` to map
  the 6 Voigt components xx,yy,zz,xy,xz,yz.

## Diseño / decisiones

- `db_version_copy(VERSION_NORMAL, VERSION_PRINT)` +
  `renumbering(VERSION_PRINT, NO, 0, 0, ...)` before reading `NODE_DOF`
  (same as print_vtk).
- The coordinates file `coord.*` is written only the first time (guarded
  by `std::ifstream` existence check), so it does not grow per time step.
- File naming follows the switch like `print_frd`/`print_gmsh`:
  `dof.<index>` for `-yes`/`-separate_index`, `dof.<n>` for
  `-separate_sequential` (static counter).

## Detalles

- Dof detection uses `dof_label` + `dof_scal_vec_mat`
  (`-SCALAR`/`-VECTOR`/`-MATRIX`), same as `print_vt` / `print_tb`.
- The `dof.<index>` file grows as a time series (append per call), one
  block of all dofs per time step. The time value is NOT written (the
  caller correlates by block position, or uses the SQLite/CSV outputs).

## Pendiente

- `control_print_dof_id` is implemented (node number column; default
  `-yes`, see [control_print_dof_id](control_print_dof_id.md)).
  `_smooth_dof`, `_smooth_n` and `_line` variants are not implemented.
- No time column: reconstructing the time of a block requires knowing the
  step sequence.
- `-separate_sequential` numbering uses a static counter (resets per run).
