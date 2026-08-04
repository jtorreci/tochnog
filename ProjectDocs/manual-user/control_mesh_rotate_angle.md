# control_mesh_rotate_angle

## Description

Rotates the 2D mesh around the z-axis by an angle given in degrees. Every node
`(x, y)` is rotated; element connectivity is not modified. Useful to orient a
mesh before solving.

NOTE: this is a PLANAR rotation. It does not generate 3D elements; the 2D to 3D
`control_mesh_rotate` remains unimplemented.

## Usage

```
control_mesh_rotate_angle <index> angle
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `angle` | Rotation angle around the z-axis, in degrees. |

## Example

Rotate the mesh 90 degrees counter-clockwise:

```
control_mesh_rotate_angle 0 90.
```
