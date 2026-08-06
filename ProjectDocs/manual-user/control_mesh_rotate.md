# control_mesh_rotate

## Description

Rotates a flat 2D mesh to 3D: the mesh is swept around the original y-axis,
which becomes the new z direction. Each `-tria3` becomes a `-prism6` and each
`-quad4` becomes a `-hex8`, so the model ends with the same number of elements
but every node is duplicated (original plus rotated copy). Useful to build a
volume of revolution from a 2D profile without remeshing.

## Usage

```
control_mesh_rotate <index> n
```

`n` is the number of elements in the rotational direction over 360 degrees.
Currently only `n = 1` is implemented, which creates a single rotational
segment.

The model must run in 3D (`number_of_space_dimensions 3`), and the integration
points must be enough for the resulting 3D elements (8 for `-hex8`, 6 for
`-prism6`).

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `n` | Number of elements in the rotational direction over 360°. Only `1` implemented. |

## Example

Revolve a 2D profile one rotational segment around the y-axis:

```
number_of_space_dimensions 3
number_of_integration_points 8
end_initia

node 1  1. 0. 0.
node 2  2. 0. 0.
node 3  2. 1. 0.
node 4  1. 1. 0.
element 1  -quad4 1 2 3 4

control_mesh_rotate 0 1
end_data
```

After this record the `-quad4` becomes a `-hex8`: the original 4 nodes plus 4
rotated copies (`x,y,z -> z,y,-x`), and the 2D source element is deleted.
