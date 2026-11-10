# prism12

## Description

`-prism12` is the 12-node prism (wedge) of the Professional manual (6.413:
"-prism12 (12 noded prismatic)"). It is quadratic in the triangle and
linear along the prism axis: the product of the 6-node triangle `-tria6`
and a 2-node line, that is a tria6 at the bottom and a tria6 at the top,
without a middle layer. Use it where the field varies quadratically in
the cross-section and slowly along the axis (extruded 2D meshes, layered
soils); the [prism18](prism18.md) adds the middle layer and is quadratic
along the axis too.

`control_mesh_extrude` turns a `-tria6` mesh into `-prism12` elements with
`control_mesh_extrude_element -prism12` (see
[control_mesh_extrude_element](control_mesh_extrude_element.md)).

## Node ordering

Two tria6 layers:

1. nodes 1-6: the bottom triangle in the `-tria6` order (corner, mid,
   corner, mid, mid, corner: corners 1, 3, 6; node 2 between 1 and 3,
   node 4 between 1 and 6, node 5 between 3 and 6),
2. nodes 7-12: the top triangle, the same order (node 6+k above node k).

For a unit prism with base corners (0,0), (1,0), (0,1) and height 1:

```
node  1 0   0   0      node  7 0   0   1
node  2 0.5 0   0      node  8 0.5 0   1
node  3 1   0   0      node  9 1   0   1
node  4 0   0.5 0      node 10 0   0.5 1
node  5 0.5 0.5 0      node 11 0.5 0.5 1
node  6 0   1   0      node 12 0   1   1
element 1 -prism12 1 2 3 4 5 6 7 8 9 10 11 12
```

The Professional 02-08-2026 uses the same order (its `element_intpnt_h`
equals the tria6 x line product to 6e-14) and so does its
`control_mesh_extrude_element -prism12`.

## Usage

Place the elements in the data part like any other element
(`number_of_space_dimensions 3`). The element works in the mixed
u-sigma formulation (`materi_velocity` / `materi_displacement` /
`materi_stress`) with `post_point`, `post_node`, `post_element_force`,
boundary conditions on nodes and geometries, gravity, etc. It is a
continuum element: an element group with `group_interface` is refused
(see Limitations).

## Integration and nodal stresses

- 14 integration points: the 7-point tria6 rule (manual 6.623: "For -tria6
  elements a seven-point scheme will be used for -maximal") x 2 Gauss
  levels along the axis (6.623: "the number of integration points in a
  direction is set equal to the number of nodes in the direction"), the
  triangle points fastest; method 0. The rule integrates the volume of any
  prism12 exactly.
- `group_integration_points -minimal` keeps the 14 points (see the
  differences below); `group_integration_method` has no effect (as for the
  prism18; the Professional also ignores `-lobatto`).
- Nodal stresses: per element the least-squares plane of the 14 point
  stresses in physical coordinates, arithmetic mean over the elements (as
  the prism15 and prism18).

## Differences with the Professional version

- Default and `-maximal` rule: the Professional 02-08-2026 stops on a
  `-prism12` with "Error detected for data item : element_intpnt_iso_coord"
  (also with `-lobatto` or `-reduced`); it runs only with `-minimal`. The
  GNU uses the 14 points of manual 6.623 (DIVERGENCES PRISM12-DEFAULT,
  class D).
- `group_integration_points -minimal`: the GNU keeps the 14 points. The
  Professional uses 8 points (the 4-point tria6 rule with the -27/48
  centroid weight x 2 Gauss levels); that rule gives an indefinite element
  stiffness (3 negative eigenvalues on a regular prism of height 0.2, up to
  4 on distorted elements) and the Professional stops on a single
  fixed-base prism of height 0.2 ("Error with pardiso solver"). Same
  decision as the prism15 (C82) and the prism18 (C132): DIVERGENCES
  PRISM12-MINIMAL, class A.
- `post_point` (and every interpolation through `point_el`): the GNU
  interpolates with the 12 shape functions; the Professional does not
  (u = 0.01 x^2 (1+z) at (0.25, 0.25, 0.5): GNU 9.375e-4 = analytic,
  Professional 5e-3). DIVERGENCES PRISM12-POST-POINT, class D.
- Interface: manual 6.624 also lists `-prism12` for `group_interface`; the
  GNU refuses it (DIVERGENCES PRISM12-INTERFACE, class B).

## Limitations

- 3D only.
- No interface `-prism12`: an element in a `group_interface` group stops
  with "Error: an interface -prism12 is not available yet; element ...".
  `control_mesh_convert` does not create interface prism12 elements from
  `-tria6` or `-prism15` either.
- Like the prism15 and prism18, no face loads through `area()`
  (`force_edge` and the edge-normal records skip the element), and the
  mesh tools that list element types (refinement, split, mpc polishing,
  geometry_moving sampling) do not handle it (BACKLOG C133).

## Tests

`validation-suite/test-2026/`: `prism12_gauss1` (unit prism, self-weight:
the 14 points, reaction -rho g V), `prism12_minimal1` (`-minimal` and
`-lobatto` give the same results), `prism12_volume1` (element with a
curved edge and a sloped top, analytic volume 2.3002 through the
self-weight reaction), `prism12_patch1` (linear patch test, 8 elements in
two layers, 5 free nodes), `prism12_post_point1` (a field quadratic in
x, y and linear in z: post_point and the stresses at the 14 points exact),
`prism12_post_element_force1` (section statics -123.4 / 1000 / -100000),
`prism12_interface_error1`.

## Related

- [prism18](prism18.md) — the 18-node wedge (tria6 x 3-node line).
- [prism15](prism15.md) — the 15-node wedge.
- [element_3d](element_3d.md) — the volume elements.
- [control_mesh_extrude_element](control_mesh_extrude_element.md).
