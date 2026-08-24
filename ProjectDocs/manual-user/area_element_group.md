# area_element_group family

## Description

Assign an element group to elements based on their location
(manual Professional 6.1-6.6). The record layout is
`-geometry_entity index element_group`; the elements matching the
selection get their `element_group` record rewritten.

Companions (same index):

- `area_element_group_method` — `-all` (default; all element nodes
  inside the geometry), `-any` (any node inside),
  `-any_but_not_all`, or a positive integer N (at least N nodes
  inside).
- `area_element_group_element` — only elements with this element name
  (e.g. `-quad4`); default all elements.
- `area_element_group_node` — select the elements by DIRECT global node
  numbers instead of a geometry (an element matches when the method
  criterion applies to its membership in that node list).
- `area_element_group_interface` — `-no` (the only legal value in
  Professional): interface elements are excluded from the re-grouping.
  The GNU also accepts `-yes` (include them), useful with the Carril A
  interface elements.
- `area_element_group_time` — `-yes`: re-evaluate the records at all
  times (default: only when the mesh changes).

The `_sequence` subfamily (manual 6.7-6.13) changes the element group
IN TIME: with `area_element_group_sequence_time` +
`area_element_group_sequence_element_group` (the Professional name; the
legacy GNU keyword `area_element_group_sequence_elementgroup` without
underscores is equivalent) group_i becomes active at time_i. Selection
by `area_element_group_sequence_geometry` (+ `_geometry_method`,
same methods as above) and/or `area_element_group_sequence` (element
numbers), with `area_element_group_sequence_element` filtering by
element name and `area_element_group_sequence_interface` excluding
interface elements (`-no`, as above).

Typical use: staged excavation/construction — switch material models
zone by zone over time (see also `control_mesh_change_element_group`;
use `control_reset_geometry` to reset stresses at the switch).

## Example

```
geometry_brick 1  0.5 0.5 0.  1.5 1.5 3.0  0.01

area_element_group        1  -geometry_brick 1  1
area_element_group_method 1  -any
area_element_group_time   1  -yes
```

```
area_element_group_sequence_geometry         0  -geometry_brick 1
area_element_group_sequence_geometry_method  0  -any
area_element_group_sequence_time             0  0.  0.1
area_element_group_sequence_element_group    0  1  2
```

Test `aeg_node`: two columns, the left one (any node inside the brick)
switches to a group with doubled Young's modulus — its sigyy doubles
(-4.0 vs -2.0). Test `aeg_seq`: the switch happens at t=0.1; with three
steps the final stress mixes both moduli exactly (-2.0).
