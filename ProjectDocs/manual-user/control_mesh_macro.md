# control_mesh_macro (node and element numbering)

## Description

`control_mesh_macro index -bar|-rectangle|-brick n_x [n_y [n_z]]` with
`control_mesh_macro_parameters` generates a structured mesh of `-bar2`,
`-quad4` or `-hex8` elements. This page documents the numbering of the
generated nodes and elements, which matters whenever an input refers to
node or element numbers (`target_item -node_dof n`, node ranges, prints).

Since F4 T8 (2026-10-10) the numbering is the one of the Professional
02-08-2026: x varies fastest, then y, then z, for the nodes, for the
elements and for the local nodes of each element (tensor order).

## Usage

```
control_mesh_macro            <index>  -rectangle <group> n_x n_y
control_mesh_macro_parameters <index>  x_c y_c l_x l_y
```

`n_x`, `n_y` are node counts per direction; the block is centred at
`x_c y_c` with lengths `l_x l_y`.

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `n_x n_y [n_z]` | Number of nodes in x, y, z. |
| `x_c y_c [z_c]` | Centre of the block. |
| `l_x l_y [l_z]` | Lengths of the block. |

## Example

```
control_mesh_macro 0   -rectangle 0 3 3
control_mesh_macro_parameters 0   1. 1. 2. 2.
```

Nodes 1, 2, 3 lie at y = 0 (x = 0, 1, 2), nodes 4, 5, 6 at y = 1; element
1 is `-quad4 1 2 4 5`, element 2 `-quad4 2 3 5 6`.

## Differences with the Professional version

None left for `-bar`, `-rectangle`, `-brick` (verified with probes). Before
F4 T8 the GNU numbered x slowest (node 2 of the example at 0 1); the
Professional `mesh_activate1` target on `node_dof 5` then read another
node. Tests: `validation-suite/test-2026/mesh_macro_numbering1` (2D),
`mesh_macro_numbering2` (3D).
