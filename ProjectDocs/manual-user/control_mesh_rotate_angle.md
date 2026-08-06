# control_mesh_rotate_angle

## Description

Rotates the flat 2D mesh around the z-axis by an angle given in degrees. Every
node `(x, y)` is rotated and only the node coordinates are modified; element
connectivity is not touched. Useful to orient a mesh before solving.

NOTE: this is a PLANAR rotation. It moves nodes in the plane; it does NOT
generate 3D elements. The 2D to 3D sweep is done by `control_mesh_rotate`.

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
