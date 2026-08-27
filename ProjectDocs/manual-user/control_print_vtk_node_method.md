# control_print_vtk_node_method

## Description

`control_print_vtk_node_method` selects which node coordinates are
written to the `POINTS` block of the `.vtk` files of
`control_print_vtk` (same `icontrol` index):

- `-node`: the stored node coordinates (the `node` records).
- `-node_start_refined`: the `node_start_refined` coordinates when they
  exist (initial / refined mesh), otherwise the stored coordinates.
  **Default** (per the gid family manual 6.307).
- `-node_deformed_mesh`: the deformed coordinates (stored coordinates +
  the nodal displacement). Only meaningful when `materi_displacement`
  is initialized; without it the stored coordinates are written.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_print_vtk` record:

```
control_print_vtk               41  -yes
control_print_vtk_node_method   41  -node_deformed_mesh
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `41`      | Index of the control record. Must match a `control_print_vtk` index. |
| node_type | `-node`, `-node_start_refined` (default), `-node_deformed_mesh`. |

## Example

A Total Lagrange column compressed with `vely = -0.01` (displacement
`disy = -0.002` after 0.2 time units):

```
control_print_vtk               40  -yes
control_print_vtk_node_method   40  -node
control_print_vtk               41  -yes
control_print_vtk_node_method   41  -node_deformed_mesh
```

`tn40.vtk` writes the top node at `y = 1` (stored); `tn41.vtk` writes it
at `y = 0.998` (deformed).
