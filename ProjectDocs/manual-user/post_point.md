# post_point

## What it does

`post_point index x y [z]` asks for the solution at an arbitrary point of
the mesh. At every printed step the GNU finds the element that contains
the point and interpolates the nodal dofs (`node_dof`) with the shape
functions of that element at the point. The result is the record
`post_point_dof index ...`, in the order of `dof_label`.

## Why it is useful

It gives results at a point that is not a node (a gauge, a sample point
of a test), and it can be used in `target_item` (`-post_point_dof index
-dof`) to check an input.

## Syntax

```
post_point  index  x y [z]
```

- `index`: the post point number.
- `x y [z]`: the point coordinates, one value per space dimension. They
  are matched against the initial (start) node coordinates of the mesh.

See also [post_point_move](post_point_move.md) (a post point that
follows the material).

## Example

```
post_point   1  1.3 0.37 0.61
target_item  1  -post_point_dof 1 -disy
target_value 1  -1.805292272921e-05 1.e-11
```

(`validation-suite/test-2026/post_point_extrude_hex27_1.dat`.)

## GiD plots

With `control_print_gid` the post points are drawn in the GiD mesh, as the
Professional 02-08-2026 does: a Point mesh `PostPoint` with one node per
active post point at its coordinates, numbered after the highest node, and
one Point element per node, numbered after the highest element, with the
material number 1000 plus the post point index. Every post point is
written, also one outside the mesh. The post point nodes carry no nodal
results. With slip surfaces (`control_safety_slip`) the post point nodes and
elements follow the slip surface nodes and elements.
(`control_print_gid_post_point1`, `control_print_gid_post_point2`.)

## Differences with the Professional version

- GiD mesh with slip surfaces (POST-POINT-GID-ELEMENT-NUMBERS, deliberate
  divergence, project owner 2026-11-25): the Professional numbers the post point elements
  after the solid mesh, on the numbers of the first slip surface elements
  (duplicate element numbers); the GNU numbers them after the slip
  elements. Without slip surfaces both write the same mesh.

- Curved or distorted quad9 / hex27 elements (deliberate divergence C78,
  project owner 2026-10-13): the GNU locates the point with the full
  isoparametric map of the element, so a post point on a node returns
  the nodal value. The Professional 02-08-2026 locates it with the
  corner nodes only (bilinear / trilinear map) and then interpolates with
  the quadratic shape functions: on a curved quad9, a post point on the
  centre node returns velx 9.437e-7 for a nodal 1.898e-7. Straight-edged
  elements with centred mid nodes give the same result in both.
- `tet10`: the GNU interpolates quadratically, the Professional linearly
  from the corners (deliberate divergence TET10-PPOINT).
- quad9 -> hex27 extruded meshes: equal to the Professional since F4 T6f
  C57 (2026-10-13); before, the inner-plane nodes had wrong start
  coordinates (see [control_mesh_extrude](control_mesh_extrude.md)).
