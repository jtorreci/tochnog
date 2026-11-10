# control_mesh_convert_quad9_quad6

## Description

Converts quadratic quadrilaterals to the 6-noded quadrilateral (Professional
manual 6.160, 2D calculations only): each `-quad9` loses the nodes of the
middle of the isoparametric direction `dir`, which becomes the linear
direction of the resulting `-quad6`. The nodes no element uses any more are
deleted; a node another element still uses is kept. Useful to obtain a mesh
that is quadratic along one direction and linear along the other.

## Usage

```
control_mesh_convert_quad9_quad6 index dir
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control index; the conversion runs when that index is processed (also without `control_timestep`). |
| `dir` | `-x` deletes the middle column of the quad9 (nodes 2, 5, 8 of the input order) and makes xi the linear direction; `-y` deletes the middle row (nodes 4, 5, 6) and makes eta the linear direction. Any other value is an input error. |

In a 1D or 3D calculation the record has no effect (no quad9 element exists).

## Example

```
number_of_space_dimensions 2
materi_velocity
materi_displacement
materi_stress
end_initia
node 1 0 0
node 2 1 0
node 3 2 0
node 4 0 1
node 5 1 1
node 6 2 1
node 7 0 2
node 8 1 2
node 9 2 2
element 1 -quad9 1 2 3 4 5 6 7 8 9
control_mesh_convert_quad9_quad6 5 -y
...
end_data
```

The result is `element 1 -quad6 1 2 3 7 8 9`; with `-x` it is
`element 1 -quad6 1 4 7 3 6 9`. Tests:
`validation-suite/test-2026/control_mesh_convert_quad9_quad6_x1.dat`,
`control_mesh_convert_quad9_quad6_y1.dat`,
`control_mesh_convert_quad9_quad6_error1.dat`.

## The continuum -quad6 element

`-quad6` is the 6-noded quadrilateral with 2 sides of 3 nodes of manual
6.413: `element index -quad6 n1 n2 n3 n4 n5 n6`, with n1-n2-n3 on one side
and n4-n5-n6 on the other, quadratic along xi (n1 to n3) and linear along
eta. It is a 2D element; in a 1D or 3D calculation it is an input error.
Outside an element group with `group_interface` it is a continuum element
(an interface `-quad6` keeps its meaning).

Integration: 3 x 2 Gauss points by default; `-minimal`: 2 x 1 Gauss points
(as the Professional; the rule is rank deficient, with 3 zero-energy modes
per element, so restrain or use the default); `-lobatto`: 3 x 2 Lobatto
points at the nodes; `-reduced`: the Lobatto points of the default rule. The
Jacobian of a -x converted element is negative (the volume uses its absolute
value, as the Professional).

Supported with the continuum quad6: stress recovery, `post_point` / point
location, `force_edge` on its 3-node and 2-node sides, `post_element_force`,
VTK output (`control_print_vtk`, native cell type 30, see
[control_print_vtk](control_print_vtk.md)).

Not supported yet (the element is skipped by these features): tendons
through the element, `control_mesh_cut_geometry`, `geometry_moving`, mpc
element polishing, `post_force_edge` totals and `element_normal`
consistent normals (see the developer page).

## Differences with the Professional version

- `-quad6` in a 3D calculation: the Professional runs it, the manual
  lists it in 2D only and the GNU refuses it (DIVERGENCES QUAD6-3D).
