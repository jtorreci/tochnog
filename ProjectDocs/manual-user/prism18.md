# prism18

## Description

`-prism18` is the 18-node quadratic prism (wedge) of the Professional
manual (6.413: "-prism18 (18 noded prismatic)"). It is the product of the
6-node triangle `-tria6` and a 3-node line along the prism axis: three
tria6 layers (bottom, middle, top). Compared with the 15-node
[prism15](prism15.md) it also has a node at the centre of each of the
three rectangular side faces (the mid-edge nodes of the middle tria6
layer).

`control_mesh_extrude` turns a `-tria6` mesh into `-prism18` elements
(see [control_mesh_extrude_element](control_mesh_extrude_element.md)).

## Node ordering

The manual side table (9.16: "For a prism18 element the sides have the
nodes numbers 0,1,2,3,4,5 and 12,13,14,15,16,17 and 0,1,2,6,7,8,12,13,14
and ...") and the element figure give three tria6 layers:

1. nodes 1-6: the bottom triangle in the `-tria6` order (corner, mid,
   corner, mid, mid, corner: corners 1, 3, 6; node 2 between 1 and 3,
   node 4 between 1 and 6, node 5 between 3 and 6),
2. nodes 7-12: the middle layer, the same order (corners 7, 9, 12 are the
   mid nodes of the vertical edges; 8, 10, 11 the side-face centres),
3. nodes 13-18: the top triangle, the same order.

For a unit prism with base corners (0,0), (1,0), (0,1) and height 1:

```
node  1 0   0   0      node  7 0   0   0.5    node 13 0   0   1
node  2 0.5 0   0      node  8 0.5 0   0.5    node 14 0.5 0   1
node  3 1   0   0      node  9 1   0   0.5    node 15 1   0   1
node  4 0   0.5 0      node 10 0   0.5 0.5    node 16 0   0.5 1
node  5 0.5 0.5 0      node 11 0.5 0.5 0.5    node 17 0.5 0.5 1
node  6 0   1   0      node 12 0   1   0.5    node 18 0   1   1
element 1 -prism18 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
```

The Professional 02-08-2026 uses the same order (its `element_intpnt_h`
equals the tria6 x line product to 1e-13).

## Usage

Place the elements in the data part like any other element
(`number_of_space_dimensions 3`). The element works in the mixed
u-sigma formulation (`materi_velocity` / `materi_displacement` /
`materi_stress`) with `post_point`, `post_node`, `post_element_force`,
boundary conditions on nodes and geometries, gravity, etc.

## Integration and nodal stresses

- 21 integration points: the 7-point tria6 rule (manual 6.622) x 3 Gauss
  levels along the axis, the triangle points fastest. This is the rule of
  the Professional for the default and `-maximal` (same point
  coordinates and order, method 0). `group_integration_method` has no
  effect (manual 6.620: bars, quads and hexes; the Professional also
  ignores it).
- Nodal stresses: per element the least-squares plane of the 21 point
  stresses in physical coordinates, arithmetic mean over the elements (as
  the prism15); equal to the Professional on the tests below.

## VTK output

`control_print_vtk` writes the element as the native VTK cell type 32
(VTK_BIQUADRATIC_QUADRATIC_WEDGE, 18 points: the 6 corners, the 9 mid-edge
nodes and the 3 side-face centres), so ParaView shows the quadratic
geometry.

## Differences with the Professional version

- `group_integration_points -minimal`: the GNU keeps the 21 points. The
  Professional uses 8 points (the 4-point tria6 rule with the -27/48
  centroid weight x 2 Gauss levels); that rule gives an indefinite
  element stiffness (6 negative eigenvalues on a distorted prism18) and a
  single element fails in the Professional solver ("Error with pardiso
  solver"). Same reason as the prism15 decision C82; the owner
  confirmation for the prism18 is pending (DIVERGENCES PRISM18-MINIMAL,
  BACKLOG C132).
- `post_point` (and every interpolation through `point_el`): the GNU
  interpolates with the 18 quadratic shape functions; the Professional
  interpolates linearly between the corner nodes (u = 0.01 z^2 at z 0.25:
  GNU 6.25e-4 = analytic, Professional 2.5e-3). DIVERGENCES
  PRISM18-POST-POINT, class D.
- VTK: the Professional writes the 6 corners as a linear wedge (type 13);
  the GNU writes the 18-point cell (owner rule 2026-11-09, DIVERGENCES
  VTK-PRISM18).

## Limitations

- 3D only.
- Like the prism15, no face loads through `area()` (`force_edge` and the
  edge-normal records skip the element), and the mesh tools that list
  element types (refinement, split, mpc polishing, geometry_moving
  sampling) do not handle it (BACKLOG C133).
- The 12-node `-prism12` of the manual is not available yet (F5 T4 U9c).

## Tests

`validation-suite/test-2026/`: `prism18_gauss1` (unit prism, self-weight:
points, reaction -rho g V, displacements and nodal stresses =
Professional), `prism18_volume1` (distorted element, analytic volume
1.825 through the self-weight reaction), `prism18_patch1` (linear patch
test, 4 elements, 5 free nodes), `prism18_post_point1` (quadratic field,
analytic), `prism18_minimal1`, `prism18_post_element_force1` (section
statics -123.4 / 1000 / -100000 = Professional), `prism18_vtk1` (VTK cell
type 32).

## Related

- [prism15](prism15.md) — the 15-node wedge.
- [element_3d](element_3d.md) — the volume elements.
- [control_mesh_extrude_element](control_mesh_extrude_element.md).
