# bounda_dof_cylindrical

## Description

Prescribes the velocity of a node in the radial direction relative to a line
defined by two points, `first` and `second`. The magnitude prescribed by
`bounda_dof`/`bounda_unknown` plus `bounda_time` is interpreted as a speed
perpendicular to the line axis (radial to the axis). The actual velocity
components applied to the node are

```
v_i = magnitude * (coord_i - proj_i) / r
```

where `proj` is the orthogonal projection of the node onto the line and `r` is
the distance from the node to the line. The node therefore moves radially
toward or away from the axis while keeping the prescribed radial magnitude.

In 1D only the `x` component is used; in 2D `x, y` (the line is taken
perpendicular to the plane, only its in-plane projection matters); in 3D all
three components of both points are used.

Useful for cylindrical flows around an axis (boreholes, tunnels). Requires
`materi_velocity` to be active and a `bounda_dof`/`bounda_unknown` record
prescribing the velocity on the same nodes and index.

## Usage

```
bounda_dof_cylindrical <index> <x_first> <y_first> <z_first> <x_second> <y_second> <z_second>
```

## Parameters

| Parameter     | Meaning                                                        |
|---------------|----------------------------------------------------------------|
| `index`       | Record index; matches the index of the `bounda_dof`/`bounda_unknown` record. |
| `x_first`     | x-coordinate of the first point of the line.                   |
| `y_first`     | y-coordinate of the first point of the line (2D/3D only).      |
| `z_first`     | z-coordinate of the first point of the line (3D only).         |
| `x_second`    | x-coordinate of the second point of the line.                  |
| `y_second`    | y-coordinate of the second point of the line (2D/3D only).     |
| `z_second`    | z-coordinate of the second point of the line (3D only).        |

## Example

Node 2 is forced to move radially with respect to the z-axis, defined by the
points `(0., 0., 0.)` and `(0., 0., 1.)`, with velocity magnitude given by
`bounda_time` (the `-velx` component carries the radial magnitude):

```
bounda_dof 1 2 -velx -vely
bounda_dof_cylindrical 1 0. 0. 0. 0. 0. 1.
```

For node 2 at `(1., 0., 0.)` the projection onto the z-axis is `(0., 0., 0.)`
and the applied velocity is radial to that axis.
