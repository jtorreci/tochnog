# element_intpnt_h / element_intpnt_iso_coord / element_intpnt_method / element_intpnt_npoint

## Description

Four print records (manual Professional 6.438, 6.439, 6.442, 6.443) that
show, per element, the space integration rule the calculation actually
uses:

| record | values | meaning |
|---|---|---|
| `element_intpnt_method index method` | 1 | the integration method: `-gauss`, `-lobatto` or `-reduced` for the bar, quadrilateral and hexahedron elements and the interfaces with a quadrilateral face; `0` (an integer, not a name) for the triangles, tetrahedra and prisms |
| `element_intpnt_npoint index npoint` | 1 | the number of integration points |
| `element_intpnt_iso_coord index ...` | npoint x ndim | per point the isoparametric coordinates |
| `element_intpnt_h index ...` | npoint x nnol | per point the value of the interpolation polynomial of every node of the element |

They are useful to check which rule an element gets from
`group_integration_method` / `group_integration_points` (and from the
automatic switches, for example the Gauss rule of an axisymmetric Lobatto
element touching the axis), and to interpolate nodal values to the
integration points of a post-processing script.

The records are written in every run, for every element, in the database
(`.dbs`) at the end of every time step; there is no switch. They are
print-only: they are ignored when read from an input file.

## Layout

- Point-major: first all the values of point 1, then point 2, ...
- `element_intpnt_iso_coord`: the bar, quadrilateral and hexahedron points
  are tensor points in [-1,1], xi fastest. Triangles give the area
  coordinates (L1, L2), tetrahedra (L1, L2, L3), prisms (L1, L2, zeta) with
  zeta in [-1,1].
- `element_intpnt_h`: per point the `nnol` nodal values in the node order
  of the element.
- Interface elements (`group_interface -yes`): one point per facing node
  pair (node i of side 1 and node i of side 2); `h` is 1/2 on both nodes of
  the pair. 2D sides: Lobatto points (xi, 0), method `-lobatto`. 3D
  quadrilateral faces (hex8, hex18): the Lobatto points (xi, eta, 0), xi
  fastest, method `-lobatto`. 3D triangular face (prism6): the corners
  (1,0,0), (0,1,0), (0,0,0), method `0`.
- Truss, beam, spring and contact spring elements get no records.

## Example

`element_intpnt_quad4_1.dat` (validation-suite/test-2026), one unit quad4
with the default rule, gives

```
element_intpnt_h  1  0.622008 0.166667 0.166667 0.0446582 0.166667 0.622008 ...
element_intpnt_iso_coord  1  -0.57735 -0.57735 0.57735 -0.57735 -0.57735 0.57735 0.57735 0.57735
element_intpnt_method  1  -gauss
element_intpnt_npoint  1  4
```

Other tests: `element_intpnt_tria3_1`, `element_intpnt_interface_quad4_1`,
`element_intpnt_1d1`, `element_intpnt_2d1`, `element_intpnt_3d1`,
`element_intpnt_interface2d1`, `element_intpnt_interface3d1`.

## Differences with the Professional version

- `element_intpnt_h` is not written when npoint x nnol exceeds 190
  values (the hex27 3x3x3 rule, the prism15 21-point rule, the quad16
  4x4 rule); the Professional writes it (729 values for a hex27)
  (DIVERGENCES INTPNT-H-LENGTH, class B). The other three records are
  written.
- The records show the GNU rule, so the known rule divergences appear in
  them: the prism6 interior triangle points (C75, the Professional uses
  the triangle corners) and the default `-reduced` factor 10 (C80, the
  Professional applies 1: points at the nodes).
- The records are filled at the start of every time step, after the mesh
  operations of the step; an element created at the end of the last step
  (refinement) has no records yet.
