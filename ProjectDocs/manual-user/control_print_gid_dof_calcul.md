# control_print_gid_dof_calcul

Manual Professional 6.300.

## Description

Only the listed post fields (see `post_calcul`) are written to the GiD results; `-none` writes
none. Without the record every post field is written. A name selects a post field when it is its
unknown (`-materi_stress` selects every operator of the stress), its operator (`-mises`) or a part
of its label, as in [control_print_vtk_dof_calcul](control_print_vtk_dof_calcul.md).

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_dof_calcul  index  name_0 name_1 ...
control_print_gid_dof_calcul  index  -none
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `name_i` | Unknown, operator or label of a post field, or `-none`. |

## Example

```
post_calcul -materi_strain_total -average -materi_stress -mises
control_print_gid              10  -yes
control_print_gid_dof_calcul   10  -materi_stress    ( only the mises field )
```

## Differences with the Professional version

- Professional-suspect (owner decision L8): the Professional 02-08-2026 removes a named post
  field, as if `-none` were given. The GNU follows the manual: only the named fields are
  written, as `control_print_gid_dof`. Guard test: `control_print_gid_dof_calcul1`.
- `-none` behaves as in the Professional.

## Tests

`validation-suite/test-2026/`: `control_print_gid_dof_calcul1`, `control_print_gid_dof_calcul2` (checked by `scripts/build_safe.sh` in a scratch copy).
