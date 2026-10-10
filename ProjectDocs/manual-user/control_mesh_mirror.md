# control_mesh_mirror

## Description

Mirrors the mesh about the plane x=0, y=0 or z=0. Every node is duplicated with
the chosen coordinate negated, and every element is duplicated with the new
nodes, so the model ends with the double of nodes and elements. Useful to
generate symmetric meshes from a half model without remeshing.

## Usage

```
control_mesh_mirror <index> -x|-y|-z
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `-x` | Mirror about the plane x=0. |
| `-y` | Mirror about the plane y=0. |
| `-z` | Mirror about the plane z=0. |

## Example

Mirror a 2D half model about x=0:

```
number_of_space_dimensions 2
end_initia

node 1  0.  0.
node 2  1.  0.
node 3  1.  1.
element 1  -bar2 1 2

control_mesh_mirror 0  -x
end_data
```

After this record the model has 6 nodes and 2 elements, symmetric about x=0.

## Once per control index (F4 T7)

Applied once per control index, on its first step, also when the same
index carries a multi-step `control_timestep` (as the Professional). See
[control_mesh_once_per_control](control_mesh_once_per_control.md).

## After a time step (F4 T8, BACKLOG C45)

The mirror can be applied on any control index, also after earlier time
steps. The new nodes carry the nodal dofs of their source nodes and the
new elements the integration-point state (stresses, strains, history) of
their source elements, as in the Professional 02-08-2026 (a mirrored
loaded block gives the Professional values to round-off). New nodes are numbered
`old + max_node + 1` and new elements `old + max_element + 1` (the
Professional numbering). Before F4 T8 the GNU stopped with a
`node_dof_start_refined` error.

Tests: `validation-suite/test-2026/mesh_copy_timestep1`,
`mesh_mirror_timestep1` (targets = Professional values).

Nodal dofs are copied as they are (no sign change of the mirrored
velocity components) and nodes on the mirror plane are not merged, as in
the Professional.
