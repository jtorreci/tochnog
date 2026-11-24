# control_print_node

## Implementación

- New file `print_node.cc`: `print_node( icontrol, ival, nval )`.
  Invoked from the control loop in `top.cc`, INSIDE the
  `control_print_frequency_allowed` gate (same as every other
  `control_print_*`):
  ```
  if ( frequency_allowed && db_active_index( CONTROL_PRINT_NODE, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_NODE, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
    print_node( icontrol, ival, ldum );
  }
  ```
  `ival[0]` is the data item (negative record number, e.g. -NODE_DOF);
  `ival[1..]` are the selected parts.
- Registered in `database.cc`: `CONTROL_PRINT_NODE` (INTEGER, variable
  length). Enums `CONTROL_PRINT_NODE*` (6) in `tochnog.h` /
  `tochnog-mod.h` in sync. `print_node.o` added to the makefile.
- Validation: the data item must be a record whose `db_name` starts
  with `node` and of type INTEGER or DOUBLE_PRECISION; otherwise
  `db_error(CONTROL_PRINT_NODE)`.

## Diseño / decisiones

- **One file per selected part** (decision, manual 6.330 example):
  - dof labels -> `db_name(labs(label)).<icontrol>` (e.g. `velx.10`,
    consistent with `control_print_dof_line`); the label is resolved
    against `dof_label` with `array_member` (unmatched -> error).
  - node_dof_calcul labels -> one file per matched item
    `post_calcul_names_without_extension[icalcul].<icontrol>` (e.g.
    `avel.10`). Matching (GNU adaptation, `post_calcul_label` does not
    exist): exact match on the underlying unknown
    (`post_calcul_unknown_operat[icalcul*2+0] == part`, e.g.
    `-materi_velocity`) OR substring of the item label
    (`strstr(post_calcul_names[icalcul], db_name(part))`, e.g. `-sigyy`
    matches `asigyy`) — same pattern as `control_print_vtk_dof_calcul`.
  - numeric parts -> `<record_name>_<n>.<icontrol>` (decision, e.g.
    `node_dof_0.10`); the number is the position in the record values
    (0-based).
  - NO parts -> ALL parts (`nuknwn` for node_dof, `ncalcul` for
    node_dof_calcul, `db_data_length` otherwise) — decision, the manual
    does not state a default.
- **Line format**: `x y z <value>` (1D: `x`), one line per node,
  append mode. Nodes without an active record, and parts beyond the
  per-node record length, are omitted for that file (decision).
- **Data source**: VERSION_NORMAL (no renumbering; node order =
  ascending node id). INTEGER records are printed as integers.
- The other five records are read per icontrol: see their developer
  pages.

## Detalles

- `nval<1` -> error. `nuknwn==0` models have no dofs; the routine
  returns when `max_node<0`.
- The `-all` special value is NOT supported as a part selector
  (`-all` in the parts -> error; `-all` is only a dof label of
  `control_print_dof_smooth_dof`).

## Pendiente

- The whole `control_print_node` family is dispatched only per
  `control_timestep` step; there is no `-separate_sequential` variant
  (the manual has no switch for it).
