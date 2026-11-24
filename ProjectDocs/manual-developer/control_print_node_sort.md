# control_print_node_sort

## Implementación

- CONTROL INTEGER record (1 value: sort method), read in
  `print_node()` (`print_node.cc`). `data_required =
  CONTROL_PRINT_NODE` in `database.cc`.
- Sort method resolution: with
  [`control_print_node_angular`](control_print_node_angular.md) only
  `-angle` is accepted (`sort_axis=3`); without angular, `-x`
  (`sort_axis=0`), `-y` only for `ndim>=2` (`sort_axis=1`), `-z` only
  for `ndim==3` (`sort_axis=2`); anything else -> `db_error`.

## Diseño / decisiones

- **Collect+sort pattern** from `print_interface_stress.cc:72-89`: the
  printable nodes (after the geometry filter) are collected into
  `node_list[]` with a sort `key[]` (the coordinate of the sort axis,
  or the angle with angular), and an INDEX array is sorted ascending
  with the same insertion sort; the files are written in that order.
  The sort is stable (equal keys keep node order).
- The key is the SAME value printed in the first column (coordinate or
  angle), so sorting matches the visible column.
- Without the record the lines are written in ascending node order
  (the natural collection order).
- The manual's mention of "control_print_node_method" is a TYPO for
  `control_print_node_angular` (6.331); it is not a keyword and is not
  implemented.

## Detalles

- One sorted order serves ALL files of the call (same node set, same
  key). The zero filter is applied per file AFTER the sort (a
  zero-suppressed line simply drops out of the sorted sequence).

## Pendiente

- None.
