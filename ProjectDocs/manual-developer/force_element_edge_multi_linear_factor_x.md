# force_element_edge_multi_linear_factor_x

## Implementación

- **Logic**: applied inside `force_factor()` in `force.cc`, at the end of
  the function. After the polynomial factor of the regular
  `force_element_edge_factor` is computed, if
  `force_element_edge_multi_linear_factor_x` is active for the same force
  index, the factor is multiplied by the multilinear value evaluated at
  `coord[0]` (x).
- **Keyword**: `force_element_edge_multi_linear_factor_x` (data_class
  FORCE, type DOUBLE_PRECISION, data_length DATA_ITEM_SIZE, fixed_length 0)
  registered in `database.cc`, with `data_required = FORCE_ELEMENT_EDGE`.
- **New enum**: `FORCE_ELEMENT_EDGE_MULTI_LINEAR_FACTOR_X` in
  `tochnog.h` / `tochnog-mod.h` (kept in sync), placed right after
  `FORCE_ELEMENT_EDGE_FACTOR`.
- **Multilinear evaluation**: pairs `x_i factor_i`. For each segment
  `[x_i, x_{i+1}]` linear interpolation; the last point extends the last
  factor to +inf; outside the first-to-last range the factor is 0
  (matches the Professional behavior).

## Diseño / decisiones

- Kept inside `force_factor()` so it applies to all callers of the edge
  force path (`area.cc`, `force.cc`) without duplicating the polynomial
  factor logic.
- Guarded by `factor_name == FORCE_ELEMENT_EDGE_FACTOR`, so it only
  affects edge forces (not normal/volume).

## Detalles

- The multiplication happens in-place: `factor *= mlfactor`.
- Only the x-coordinate is used (`coord[0]`), per the keyword name.

## Pendiente

- Only the edge-force variant is implemented. The Professional family also
  has `force_edge_node_factor` and per-node factors; those are separate
  keywords not covered here.
- In 1D the coordinate is `coord[0]` as well (same path).
