# control_mesh_extrude

## Description

Extrudes a flat 2D mesh (at z=0) to 3D along the z-axis. Each 2D element is
swept through a stack of layers, generating one 3D element per 2D element per
layer: each `-tria3` becomes a `-prism6` and each `-quad4` becomes a `-hex8`.
Useful to build a solid volume from a 2D profile without remeshing.

## Usage

```
control_mesh_extrude <index> z0 z1 z2 ...
```

`z0`, `z1`, `z2`, ... give the z-coordinate of each layer boundary, in order.
With `n` values the mesh is extruded into `n` layers; one 3D element is created
per 2D element and layer, so the model ends with `n` times the original number
of elements and `n` copies of every node (one copy per layer boundary).

The model must run in 3D (`number_of_space_dimensions 3`), and the integration
points must be enough for the resulting 3D elements (8 for `-hex8`, 6 for
`-prism6`).

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `z0 z1 z2 ...` | z-coordinates of the layer boundaries; one layer is created between consecutive values. |

## Example

Extrude a single `-quad4` with one layer up to z=1:

```
number_of_space_dimensions 3
number_of_integration_points 8
end_initia

node 1  0. 0. 0.
node 2  1. 0. 0.
node 3  1. 1. 0.
node 4  0. 1. 0.
element 1  -quad4 1 2 3 4

control_mesh_extrude 0 1.
end_data
```

After this record the `-quad4` becomes a `-hex8`: the original 4 nodes plus 4
copies at z=1, and the 2D source element is deleted. For a mesh with nodes
`1..n`, the copy of node `inod` in layer `layer` is at
`inod + (max_node+1)*(layer+1)`, i.e. the copies start right after
`max_node+1`.
