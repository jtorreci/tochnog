# force_gravity_geometry

## Description

`force_gravity_geometry geometry_item_name geometry_item_index` (manual
6.499, no index): only elements inside the geometry get the
`force_gravity` body force. Without the record every element gets it.

## Effect on the calculation

- An element is inside when ALL its nodes are inside the geometry (initial
  node coordinates). An element half inside gets no gravity, even when its
  centre is inside (Professional 02-08-2026).
- It applies to continuum (`-materi`) and truss elements. Other uses of
  `force_gravity` (nodal masses, hydrostatic water pressures,
  `post_element_force`) are not restricted.
- It needs `force_gravity`.

## Example

```
force_gravity 0. -10.
geometry_quadrilateral 1  -0.5 -0.5  1.45 -0.5  -0.5 1.5  1.45 1.5  1.e-3
force_gravity_geometry -geometry_quadrilateral 1
```

Mind the `geometry_quadrilateral` point order (manual 6.544: `0 0  1 0
0 1  1 1`, as the nodes of a quad4).

Self-checking inputs in `validation-suite/test-2026`:
`force_gravity_geometry1` (one of two columns), `force_gravity_geometry2`
(none), `force_gravity_geometry3` (trusses).
