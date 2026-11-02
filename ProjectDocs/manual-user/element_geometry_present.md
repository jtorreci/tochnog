# element_geometry_present + print_element_geometry_present

## Description

`element_geometry_present` lists for every element the geometries in
which the element is present — a diagnostic record to check that the
geometries include exactly the elements you expect (manual Professional
6.427). An element is present in a geometry when ALL its nodes are in
that geometry. The filling is switched on with
`print_element_geometry_present -yes` (6.975, default `-no`).

`print_element_geometry_present_node_type` (6.976) selects the node
coordinates of the check: `-node_start_refined` (default), `-node` or
`-plus_displacement` (node plus nodal displacements). The per-geometry
`geometry_node_type` records override it for their own geometry.

The record is filled at the start of every time step, and every step
overwrites the list of the previous one (as `node_geometry_present`).

## Usage

```
print_element_geometry_present -yes
print_element_geometry_present_node_type [-node_start_refined | -node | -plus_displacement]
```

The stored record is a list of `(geometry name, geometry index)` pairs:

```
element_geometry_present 1 -geometry_quadrilateral 7
```

## Example

`element_geometry_present1.dat` (validation-suite/test-2026): a 4 x 4
quad4 mesh of unit elements and

```
print_element_geometry_present -yes
geometry_quadrilateral 7 -0.01 -0.01 1.01 -0.01 -0.01 1.01 1.01 1.01 0.001
geometry_line 8 0 0 1 0 0.01
```

gives only `element_geometry_present 1 -geometry_quadrilateral 7`: the
line holds two nodes of element 1, not all four.

`element_geometry_present_plus_displacement1.dat`: one quad4 translating
at velx 0.2 with `-plus_displacement`; at the start of the last step the
displaced element lies in the quadrilateral (x 0.05 to 1.15) and the
record is filled, while the default `-node_start_refined`
(`element_geometry_present_start1.dat`) gives no record.

## Differences with the Professional version

- Same records and same rule (all nodes present) as the Professional
  02-08-2026.
- The Professional rejects `-plus_displacement` in
  `print_element_geometry_present_node_type` ("I do not know :
  -plus_displacement"), although its manual lists it; the GNU accepts it
  as the manual says (DIVERGENCES EGP-PLUS-DISP). Before F5 T2 V3 the GNU
  rejected the word too, in every record that lists it.
