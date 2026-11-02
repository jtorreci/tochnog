# element_geometry_present + print_element_geometry_present (developer)

## Where

- `geometry.cc`, `void element_geometry_present_calculate( void )`: the
  element counterpart of `node_geometry_present_calculate()`. For every
  active element and every active geometry of the scanned entity list,
  `geometry()` is called for each element node (stops at the first node
  outside); when all are found the pair `(-entity, index)` is appended.
  The record is PUT when the list is not empty, deleted otherwise.
- `top.cc`, `step_start()`: called right after
  `node_geometry_present_calculate()`.
- `database.cc`: `PRINT_ELEMENT_GEOMETRY_PRESENT`,
  `PRINT_ELEMENT_GEOMETRY_PRESENT_NODE_TYPE` (INTEGER, length 1, no index)
  and `ELEMENT_GEOMETRY_PRESENT` (INTEGER, `DATA_ITEM_SIZE`, variable
  length, class ELEMENT).
- Enums appended before `LAST_DUMMY` in `tochnog.h` / `tochnog-mod.h`.
- `database.cc`: `name[PLUS_DISPLACEMENT] = "plus_displacement"`. The enum
  value existed (handled in `geometry()`), but without a name the parser
  rejected `-plus_displacement` in `geometry_node_type`,
  `print_node_geometry_present_node_type` and this record.

## Implementation details

- Entity list: the same 15 geometry entities as the node version, in
  alphabetical record order (geometry_list included).
- Node type: `-print_element_geometry_present_node_type` (positive enum
  value) passed to `geometry()`; default `NODE_START_REFINED`. The
  per-geometry `geometry_node_type` overrides it inside `geometry()`.
- More than `DATA_ITEM_SIZE/2` present geometries for one element stop
  with `db_error`.

## Pending refactorings

- The entity list is duplicated from `node_geometry_present_calculate()`;
  a shared static table would keep both scans in sync.

## Verification

Professional 02-08-2026 probes (4 x 4 quad4): quadrilateral around
element 1 -> `element_geometry_present 1 -geometry_quadrilateral 7`, also
with `-node`; a quadrilateral over half of element 1 or a line through
two of its nodes -> no record. Tests `element_geometry_present1`,
`element_geometry_present_node_type1`, `element_geometry_present_partial1`,
`element_geometry_present_plus_displacement1`, `element_geometry_present_start1`.
