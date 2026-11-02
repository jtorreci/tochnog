# print_vtk_calculation / print_gmsh_calculation / print_tecplot_calculation

## Description

End-of-run output switches (manual Professional 6.999, 6.991, 6.998):
with `-yes` the vtk, gmsh or tecplot file is written once more at the end
of the calculation, with the final state; with `-no` nothing is written at
the end. They do not need any `control_print_*` record.

- `print_vtk_calculation -yes`: writes `<base>.vtk` (input `excavation.dat`
  -> `excavation.vtk`), with the same content as a `control_print_vtk`
  file (the global `print_vtk_*` defaults apply).
- `print_gmsh_calculation -yes`: appends the final state to `<base>.msh`
  (the file of `control_print_gmsh -yes`; the mesh is written first when
  the file is new). When a `control_print_gmsh -yes` already printed the
  final time, the state is not repeated.
- `print_tecplot_calculation -yes`: appends one zone with the final state
  to `tecplot.plt` (the file of `control_print_tecplot`), again skipped
  when the final time is already in the file.

Default `-no` for the three switches (see the differences below).

## Usage

```
print_vtk_calculation      -yes
print_gmsh_calculation     -yes
print_tecplot_calculation  -yes
```

| Parameter | Meaning |
|---|---|
| switch | `-yes`: write the file at the end of the run. `-no` (default): do not. |

## Example

```
control_timestep 0 0.1 0.5
print_vtk_calculation -yes
```

writes `<base>.vtk` with the state at t = 0.5. Tests
`validation-suite/test-2026/print_vtk_calculation1`,
`print_gmsh_calculation1/2`, `print_tecplot_calculation1/2`,
`print_calculation_default1`.

## Differences with the Professional version

- Default: the manual gives `-yes` for `print_tecplot_calculation` and
  `print_vtk_calculation`; the Professional 02-08-2026 writes these files
  only with `-yes`. Project-owner decision 2026-11-02 ("vamos con pro, y lo
  dejamos documentado"): the GNU default is `-no`, as the Professional
  (DIVERGENCES `PRINT-CALC-DEFAULT`, class A; SEGUIMIENTO Descarte).
- File names: the Professional appends the tecplot zones to `<base>.plt`;
  the GNU uses `tecplot.plt`, the file of its `control_print_tecplot`
  (BACKLOG C99, class B). The vtk end file `<base>.vtk` has the
  Professional name.
- The layout of the files follows the GNU writers (BACKLOG C27).
