# post_force_edge_summed

## Implementación

- **Keyword**: `post_force_edge_summed` in `database.cc` (DOUBLE_PRECISION,
  `data_length = ndim`, data_class POST, `no_index = 1` - stored at index
  0, printed without index like the Professional .dbs), enum
  `POST_FORCE_EDGE_SUMMED` in `tochnog.h`/`tochnog-mod.h` between
  `POST_ERROR_RESULT` and `POST_GLOBAL` (kept in sync).
- **Cálculo**: `post_force_edge_summed_calculate()` in `post.cc`, called
  at the top of `post()` (every step_close, BEFORE the control_print
  section prints the record and before the exit_tn target check). It is
  computed ONLY on demand: the function scans the `target_item` records
  and the `control_print` item lists for the keyword and returns without
  doing anything otherwise (zero impact on the other tests; the mesh
  scan per step would otherwise slow every force_edge run).
- **Integración**: mirrors the load assembly of `area()`
  (`force_element_edge`): for every active `force_edge` record and every
  element side whose nodes ALL lie on the record's geometry entity (or
  in its node-area list) and that passes the element/group/side/node
  restrictions (`force_element_edge_element(_group/_side/_node/
  _element_node)`), the total force vector accumulates per side node
  `w_lobatto * ar * load * factor * node_factor * values[dir]` with the
  SAME quadrature that distributes the traction to the element nodes:
  - 2D: side length (first-last side node) x Lobatto weights of the side
    (quad9: 1/6, 4/6, 1/6);
  - 3D: TET4 Heron triangle x 1/3; HEX8/HEX27 face = triangle fan
    (0,1,2)+(1,2,last) x tensor Lobatto 2x2 / 3x3.
  Border tables: local copies of the `border_nodes_*` of `area.cc`.
- **load/factor**: `force_element_edge_sine` / `_time` / `_time_file`
  resolved in the same order as `area()` (default load = 1); spatial
  factor via `force_factor()` (`force_element_edge_factor` incl. the
  multi_linear_factor_x); per-node `force_element_edge_node_factor`.
- **Direction mapping**: the record values map onto the space directions
  through the principal unknowns, the same loop as the `area()`
  assembly (`dof_principal` >= 0).
- **Verificación**: corpus `elasti6.dat` - target
  `-post_force_edge_summed 0 1` = 10 ± 1e-4: GNU value (0, 10) EXACT
  (same as the Professional .dbs of elasti6: `post_force_edge_summed 0.
  1.000000000000e+01`).
- **Gotcha**: `geometry()` receives an MDIM-sized projection buffer and
  the Lobatto ISO coordinates need their own arrays (small local buffers
  overflowed in the first draft).

## Pendiente

- The rest of the Professional `post_global` *_summed family
  (`post_bounda_force_summed`, `post_node_summed`, ...) is not
  implemented; `-post_force_edge_summed` was added for the corpus
  `elasti6` check.
- `force_edge_normal/projected/water` masters are not summed (separate
  items in the Professional manual).
