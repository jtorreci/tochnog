# control_print_vtk_dof

## Description

`control_print_vtk_dof` limits the solution fields written to the `.vtk`
files of `control_print_vtk` (same `icontrol` index) to the listed ones.
This makes the VTK files smaller, which is convenient for very large
calculations.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_print_vtk` record:

```
control_print_vtk            20  -yes
control_print_vtk_dof        20  -materi_stress
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match a `control_print_vtk` index. |
| names     | Initialisation names to write, e.g. `-condif_temperature`, `-materi_velocity`, `-materi_stress`, `-materi_strain_total`. Use `-none` to write no field. |

Without this record all solution fields are written (default).

## Example

Write only the stress field (and its derived magnitudes von Mises, Tresca,
principal stresses) to the VTK files:

```
control_print_vtk            20  -yes
control_print_vtk_dof        20  -materi_stress
```

Produces a `.vtk` file containing `TENSORS materi_stress` (plus the
component scalars and the derived magnitudes `vmises`, `tresca`,
`sig1..sig3`). Use `-materi_velocity` to write only the velocity, or
`-none` to write no nodal fields at all.
