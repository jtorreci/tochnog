# bounda_dof_radial

## Description

Prescribes the velocity of a node in the radial direction relative to a fixed
point `(x, y, z)`. The magnitude prescribed by `bounda_dof`/`bounda_unknown`
plus `bounda_time` is interpreted as a radial speed: the actual velocity
components applied to the node are

```
v_i = magnitude * (coord_i - point_i) / r
```

where `r` is the distance from the node to the point. The node therefore
moves radially toward or away from the point while keeping the prescribed
radial magnitude.

In 1D only the `x` component is used, in 2D `x, y`, and in 3D `x, y, z`.

Useful for radial flows or expansion/contraction of the mesh. Requires
`materi_velocity` to be active and a `bounda_dof`/`bounda_unknown` record
prescribing the velocity on the same nodes and index.

## Usage

```
bounda_dof_radial <index> <x> <y> <z>
```

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `index`   | Record index; matches the index of the `bounda_dof`/`bounda_unknown` record. |
| `x`       | x-coordinate of the radial center point.                       |
| `y`       | y-coordinate of the radial center point (2D/3D only).          |
| `z`       | z-coordinate of the radial center point (3D only).             |

## Example

Node 2 is forced to move radially with respect to the origin `(0., 0., 0.)`,
with velocity magnitude given by `bounda_time` (`-velx` component carries the
radial magnitude, here -1 at all times):

```
bounda_dof 1 2 -velx -vely
bounda_time 1 0. -1. 100000. -1.
bounda_dof_radial 1 0. 0. 0.
```

For node 2 at `(1., 0.)` the applied velocity is `(-1, 0)`; for node 2 at
`(0., 1.)` it is `(0, -1)` — the node moves radially toward the origin.
