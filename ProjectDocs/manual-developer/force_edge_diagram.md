# force_edge_diagram — developer

## Files and functions

- `database.cc` — `FORCE_ELEMENT_EDGE_DIAGRAM` (`force_element_edge_diagram`,
  DOUBLE, `DATA_ITEM_SIZE`, class FORCE, requires `FORCE_ELEMENT_EDGE`). The
  Professional name `force_edge_diagram` reaches it through the
  `force_edge` -> `force_element_edge` prefix translation in `db_number`.
- `tochnog.h` / `tochnog-mod.h` — enum after `FORCE_ELEMENT_EDGE_NODE_FACTOR`.
- `area.cc` — `area()`, `FORCE_ELEMENT_EDGE` branch of the side loop.

## Implementation details

In the per-node part of the `FORCE_ELEMENT_EDGE` branch, after the node
restrictions: when the index has a diagram, `table_xy` interpolates it at
`new_dof[inol*nuknwn + dis_indx + (ndim-1)*nder]` (the current iterate of
the vertical displacement) and the result multiplies `node_factor`.
`table_xy` returns 0 outside the table. The diagram needs at least two
pairs (`db_error` otherwise) and `materi_displacement` (error message
otherwise). The record values are copied before `table_xy`, which sorts
its table in place.

Measured on the Professional 02-08-2026 (two quad4 columns E = 10 / 100,
top edge load, diagram -1 0 0 1):

- first iteration factor 1, second iteration factor f(d) of the first:
  one step gives -0.09 / -0.017851239669 / -0.0099 at the top nodes, the
  GNU values to 12 digits;
- 30 steps converge to -0.1/1.1 and -0.01/1.01 at the column ends: the
  factor is taken per node (per integration point would mix the two
  columns);
- a node below the first diagram point gets factor 0.

## History

- F5 T5 U5 (2026-11-18): added.

## Hardcoded parameters / pending refactorings

- The vertical direction is the last space dimension.
- Only `force_edge` uses the diagram; the Professional manual names no
  diagram for the other edge families.
