# print_gmsh_node_type

## Description

`print_gmsh_node_type node_type` (manual Professional 6.993): default of
`control_print_gmsh_node_type` (GNU record `control_print_gmsh_node_method`,
see [control_print_gmsh](control_print_gmsh.md)) for every gmsh file:
`-node` (default), `-node_start_refined` or `-node_deformed_mesh`
coordinates in `$Nodes`. A control record overrules it for its index. The
older name `print_gmsh_node_method` is accepted.

## Usage

```
print_gmsh_node_type -node_deformed_mesh
```

| Parameter | Meaning |
|---|---|
| node_type | `-node` (default), `-node_start_refined`, `-node_deformed_mesh`. |

Test `validation-suite/test-2026/print_gmsh_node_type1` (mesh written at the
first print, node 2 x 0.999914285714 = Professional 02-08-2026).

## Differences with the Professional version

None for the switch.
