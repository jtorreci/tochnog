# control_print_node_geometry

## Implementación

- CONTROL INTEGER record (2 values: geometry item name + index), read
  in `print_node()` (`print_node.cc`). `data_required =
  CONTROL_PRINT_NODE` in `database.cc`; `fixed_length=1` with
  `data_length=2` (the record has exactly 2 values).
- Validation: `geom_entity[0] < 0` and
  `db_data_class(geom_entity[0]) == GEOMETRY`, else `db_error`.

## Diseño / decisiones

- Filter pattern from `filter.cc:74-80` (the PRINT_FILTER geometry
  case): `geometry( inod, ddum, geom_entity, in_geometry, factor,
  normal, penetration, projection, node_type, PROJECT_EXACT,
  VERSION_NORMAL )` with `geom_entity = [name, index]` (the record
  itself).
- `node_type` = `NODE_START_REFINED` when any `NODE_START_REFINED`
  record exists (`db_max_index >= 0`), else `NODE` — same fallback as
  `control_print_dof_line` (decision; the manual does not specify the
  frame).
- The filter runs when building the printable node list (before the
  sort), so the sort order and the zero filter operate only on the
  nodes in the geometry.

## Detalles

- Nodes outside the geometry are omitted from ALL files of the call.

## Pendiente

- None.
