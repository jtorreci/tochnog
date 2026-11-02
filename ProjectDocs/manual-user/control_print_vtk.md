# control_print_vtk

## Description

`control_print_vtk index switch` (manual Professional 6.339) writes the
results in the Visual Toolkit unstructured grid format (for ParaView) at
the time steps of `control_timestep index`. The switch selects the file:

| switch | file (input `excavation.dat`) |
|---|---|
| `-yes` | `excavation.vtk`, overwritten by every print (it holds the last one) |
| `-separate_index` | `excavation_<index>.vtk` with the index in 3 digits: `excavation_000.vtk`, `excavation_012.vtk`, `excavation_100.vtk` |
| `-separate_sequential` | one file per print, with a counter in 4 digits: `excavation_0000.vtk`, `excavation_0001.vtk`, ... |
| `-no` | nothing |

The content is set by the `control_print_vtk_*` records of the same index
and the global `print_vtk_*` defaults.

## Usage

```
control_timestep   0  0.1 0.2
control_timestep  12  0.1 0.3
control_print_vtk  0  -separate_index   # excavation_000.vtk
control_print_vtk 12  -separate_index   # excavation_012.vtk
```

Tests `validation-suite/test-2026/control_print_vtk_names1`,
`control_print_vtk_separate_index1`, `control_print_vtk_separate_sequential1`.

## Differences with the Professional version

None for the names (= Professional 02-08-2026). They differ from the
manual 6.339, which writes `excavation100.vtk` and `excavation0.vtk`
(no underscore, no zero padding): owner decision 2026-11-02, DIVERGENCES
`VTK-TECPLOT-NAMES` (class A). Before F5 the GNU wrote `tn<index>.vtk` for
`-yes` and ignored the two `-separate_*` switches.
