# control_mesh_copy

## Description

Copies the mesh displaced: every node is duplicated and moved over
`move_x`, `move_y`, `move_z`, and every element is duplicated with the new
nodes, so the model ends with the double of nodes and elements. Useful to
build repetitive meshes (for example periodic patterns) without remeshing.

## Usage

```
control_mesh_copy <index> move_x [move_y [move_z]]
```

Only the displacements for the number of space dimensions of the model must be
given: `move_x` in 1D, `move_x move_y` in 2D, all three in 3D.

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `move_x` | Displacement applied to the x-coordinate of every copied node. |
| `move_y` | Displacement applied to the y-coordinate of every copied node. |
| `move_z` | Displacement applied to the z-coordinate of every copied node. |

## Example

Copy a 1D mesh displaced by 2 in x:

```
number_of_space_dimensions 1
end_initia

node 1  0.
node 2  1.
element 1  -bar2 1 2

control_mesh_copy 0  2.
end_data
```

After this record the model has 4 nodes (the originals plus copies at 2. and 3.)
and 2 elements.

## Once per control index (F4 T7)

Applied once per control index, on its first step, also when the same
index carries a multi-step `control_timestep` (as the Professional). See
[control_mesh_once_per_control](control_mesh_once_per_control.md).

## After a time step (F4 T8, BACKLOG C45)

The copy can be applied on any control index, also after earlier time
steps. The new nodes carry the nodal dofs of their source nodes and the
new elements the integration-point state (stresses, strains, history) of
their source elements, as in the Professional 02-08-2026: a copy of a
loaded block behaves exactly as the original. New nodes are numbered
`old + max_node + 1` and new elements `old + max_element + 1` (the
Professional numbering). Before F4 T8 the GNU stopped with a
`node_dof_start_refined` error.

Tests: `validation-suite/test-2026/mesh_copy_timestep1`,
`mesh_mirror_timestep1` (targets = Professional values).
