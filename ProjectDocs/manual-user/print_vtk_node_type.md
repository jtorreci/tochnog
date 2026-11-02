# print_vtk_node_type

## Description

`print_vtk_node_type node_type` (manual Professional 6.1002): default of
`control_print_vtk_node_type` (GNU record `control_print_vtk_node_method`,
see [control_print_vtk_node_method](control_print_vtk_node_method.md)) for
every vtk file: `-node`, `-node_start_refined` or `-node_deformed_mesh`
coordinates in `POINTS`. A control record for a control index overrules it.
The older name `print_vtk_node_method` is accepted and stored as
`print_vtk_node_type`.

## Usage

```
print_vtk_node_type -node_deformed_mesh
```

| Parameter | Meaning |
|---|---|
| node_type | `-node`, `-node_start_refined` (GNU default) or `-node_deformed_mesh`. |

Tests `validation-suite/test-2026/print_vtk_node_type1/2`,
`print_vtk_node_method1` (deformed node 2 x 0.998714285714 = Professional).

## Differences with the Professional version

None measured (default `-node_start_refined` falls back to the stored
coordinates without refinement, the Professional writes the stored ones).
