# control_print_vtk_dof_calcul

## Description

`control_print_vtk_dof_calcul` limits the POST fields (the results
computed by `post_calcul`) written to the `.vtk` files of
`control_print_vtk` (same `icontrol` index) to the listed ones. It is
the mirror of `control_print_vtk_dof` over the post_calcul block. This
makes the VTK files smaller, which is convenient for very large
calculations.

Without this record every post field is written (default). Use `-none`
to write no post field at all (the solution fields written by
`control_print_vtk_dof` are unaffected).

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_print_vtk` record:

```
control_print_vtk               34  -yes
control_print_vtk_dof_calcul    34  -materi_strain_total
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `34`      | Index of the control record. Must match a `control_print_vtk` index. |
| names     | Post fields to write. Each name is resolved like any keyword: the initialisation name of the underlying data (e.g. `-materi_stress`, `-materi_strain_total`, `-materi_velocity`, `-condif_temperature`, `-groundflow_pressure`) selects every operator computed for that data, and a component name (e.g. `-sigyy`, `-eptyy`) selects the post fields whose label contains it. Use `-none` to write no post field. |

The Professional manual refers to `post_calcul_label` for the allowed
names; that record does not exist in the GNU, where the post field
labels are the names shown in the vtk headers (e.g.
`materi_strain_total_average` for `post_calcul -materi_strain_total
-average`, label `aept`).

## Example

With

```
post_calcul -materi_strain_total -average -materi_stress -mises
control_print_vtk               34  -yes
control_print_vtk_dof_calcul    34  -materi_strain_total
```

only the `materi_strain_total_average` post field (label `aept`) is
written; the `materi_stress_mises` field (label `mises-sig`) is not.
`-none` removes both post fields.
