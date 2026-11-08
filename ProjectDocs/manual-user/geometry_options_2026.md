# Geometry options and entities of the 2026 manual

`geometry_boundary`, `geometry_circle_part`, `geometry_cylinder_part`,
`geometry_cylinder_part_start_vector`, `geometry_element_geometry`,
`geometry_element_geometry_method`, `geometry_exclude`,
`geometry_hexahedral`, `geometry_line_eps_iso`, `geometry_mpc`,
`geometry_quadrilateral_eps_iso`, `geometry_tetrahedral` (manual
Professional 6.509-6.549).

## Entities

They can be used wherever a geometry is accepted (`bounda_dof 0
-geometry_tetrahedral 1 ...`, `control_mesh_delete_geometry`, ...).

| record | values | selects |
|---|---|---|
| `geometry_circle_part index xc yc angle_start angle_end radius tolerance` | 2D | nodes within `tolerance` of the circle and with an angle (radians, from the positive x-axis, 0 to 2 pi) in the range |
| `geometry_cylinder_part index x0 y0 z0 x1 y1 z1 radius as0 ae0 [as1 ae1 ...] tolerance` | 3D | nodes within `tolerance` of the cylinder surface between the two axis points, with the angle around the axis in one of the ranges |
| `geometry_cylinder_part_start_vector index vx vy vz` | 3D | the angle is measured from this vector (its axial part is removed); only one range `0 end`, valid for angles up to `end` on both sides |
| `geometry_tetrahedral index` 4 corner points | 3D | everything inside (boundary included) |
| `geometry_hexahedral index` 8 corner points | 3D | everything inside, corners in hex8 order (`0 0 0`, `1 0 0`, `0 1 0`, `1 1 0`, then the same at the top) |

Without a start vector the cylinder axis must be along x, y or z; the
angle is measured from x towards y (axis z), from y towards z (axis x) or
from x towards z (axis y). Angles are non-negative and every end angle is
larger than its start angle.

## Options (same index as the geometry)

| record | meaning |
|---|---|
| `geometry_line_eps_iso index eps` | accept points up to `eps` (in line coordinates) beyond the ends of the line |
| `geometry_quadrilateral_eps_iso index eps` | accept points just outside the quadrilateral (tolerance of the triangle coordinates of its two triangles), like `geometry_triangle_eps_iso` |
| `geometry_boundary index -yes/-no` | only nodes on (`-yes`) / not on (`-no`) the mesh boundary (`node_boundary`) |
| `geometry_mpc index -yes/-no` | `-yes`: only nodes constrained by an mpc (the first node of an `mpc_node_number`, also the generated ones); `-no`: no restriction |
| `geometry_element_geometry index set_0 set_1 ...` | only nodes of elements with these `element_geometry` sets |
| `geometry_element_geometry_method index method` | `-any` (default): an element of the node in a listed set; `-all`: every listed set among the elements of the node; `-only`: every element of the node in a listed set |
| `geometry_exclude index -geometry_xxx gindex ...` | not in any of the listed geometries; a geometry may not exclude itself |

The node options act on node tests; a test of a point (for example a
post point) only uses `geometry_exclude`.

## Example

```
geometry_quadrilateral 10 0. 0. 1. 0. 0. 1. 1. 1. 0.01
geometry_exclude 10 -geometry_point 20
geometry_point 20 0.5 0.5 0.3
geometry_boundary 10 -no
```

Tests (validation-suite/test-2026): `geometry_eps_iso1`,
`geometry_boundary1`, `geometry_exclude1`, `geometry_exclude_error1`,
`geometry_element_geometry1`, `geometry_circle_part1`, `geometry_mpc1`,
`geometry_tetrahedral1`, `geometry_hexahedral1`, `geometry_cylinder_part1`,
`geometry_cylinder_part_error1`, `_error2`.

Errors: `Error detected for data item : geometry_exclude` (a geometry that
excludes itself, a chain deeper than 20, an odd number of values);
`Error detected for data item : geometry_cylinder_part` (no start vector
and an axis not along x, y or z; a bad angle range);
`Error, GEOMETRY_CYLINDER_PART_START_VECTOR can only be combined with one
angle range.`

## Differences with the Professional version

- On the plane through the first axis point of a `geometry_cylinder_part`
  the Professional 02-08-2026 selects other nodes (axis z from 0,0,0, range
  0-0.8: it takes 0,2,0 and 1,2,0 and drops 2,0,0); the GNU measures every
  plane as the manual says (DIVERGENCES CYLINDER-PART-START-PLANE).
- A geometry that excludes itself is an error in the GNU (manual: not
  allowed); the Professional ignores the record. A cylinder part without
  start vector and an oblique axis is an error in the GNU; the Professional
  selects no node (DIVERGENCES GEOMETRY-INPUT-CHECKS).
- `geometry_method` (6.532) is documented in
  [element_geometry_selection](element_geometry_selection.md).
