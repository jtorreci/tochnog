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

## Quadratic sections: `-quad9` to `-hex27`

A `-quad9` element is extruded into `-hex27` elements, as in the
Professional (its `test/other/force11` ring). For this branch the `z`
record gives the total extent `z0 z1` and the first value of
`control_mesh_extrude_n` gives the number of hex27 layers:

```
control_mesh_extrude   0  0. 10.
control_mesh_extrude_n 0  2 1
```

Each hex27 layer has 9 base nodes, 9 nodes on the mid-plane at the layer
mid-height and 9 top nodes. The source nodes are the base of the first
layer and the copies at `z1` the top of the last one. The inner planes
(every mid-plane and every layer top below `z1`) hold **one node per
source node, shared by all the hex27 elements that touch it**, so the
extruded mesh is connected everywhere: two hex27 layers give 5 z planes
with as many nodes each as the source mesh (Professional `force11`: 2400
nodes per plane, 12000 in all).

Until 2026-10-11 (F4 T12, BACKLOG C47) every element received its own
inner-plane nodes, so neighbouring elements were connected only at `z0`
and `z1`; `force11` came out about 8 times too soft and missed its
section-moment target (5.80223 instead of 5.76460). Now its nodal
displacements equal the Professional's and `mom1y_sig` is 5.76461.

Self-checking input: `validation-suite/test-2026/extrude_quad9_shared_nodes1.dat`
(two quad9 x two hex27 layers, a cantilever plate: 75 nodes, nodal
displacements = Professional).

## Differences with the Professional version

- Node numbering of the quad9 -> hex27 extrusion: the GNU numbers inner
  plane `q` of source node `inod` as `next + (q-1)*(max_node+1) + inod`
  (with gaps), the Professional numbers the planes compactly in z order
  (plane `q` of node `inod` is `q*n + inod`). Geometry, connectivity and
  results are the same; node-numbered output (`node_dof N`, targets on
  node numbers) differs.
- hex27 nodal stresses (`force11`): equal to the Professional (7e-7 of
  the field maximum) since F4 T6f C70 (2026-10-13, physical-coordinate
  recovery fit; before up to 3.5 % off at the loaded points). Still open:
  a `post_point` on a hex27 node (BACKLOG C57, DIVERGENCES C57).

