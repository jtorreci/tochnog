# bounda_normal

## Description

Constrains nodes to slide on a plane: they cannot move along the direction
normal to the plane, and the velocity component normal to the plane is set to
zero. The tangential components are left free.

In 3D the full normal vector `(normal_x, normal_y, normal_z)` is given; in 2D
only the two in-plane components are used; in 1D only one component is used.

Requires `materi_velocity` to be active (the projection is applied to the
velocity dofs) and a `bounda_unknown`/`bounda_dof` record on the same nodes.

## Usage

```
bounda_normal <index> <normal_x> <normal_y> <normal_z>
```

## Parameters

| Parameter     | Meaning                                                        |
|---------------|----------------------------------------------------------------|
| `index`       | Record index; matches the index of the `bounda_unknown`/`bounda_dof` record. |
| `normal_x`    | x-component of the plane normal vector.                        |
| `normal_y`    | y-component of the plane normal vector (2D/3D only).           |
| `normal_z`    | z-component of the plane normal vector (3D only).              |

## Example

Nodes may only move in the plane perpendicular to the x-axis (`velx = 0`):

```
bounda_normal 1 1. 0. 0.
```

With `bounda_dof 1 ...` (or `bounda_unknown 1 ...`) present for the same index,
the node slides freely in the y-z plane while the x-velocity is annulled.
