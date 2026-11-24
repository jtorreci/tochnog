# control_print_gid_dof

Manual Professional 6.299.

## Description

Only the listed solution fields are written to the GiD results. The names are initialisation
names (`-materi_velocity`, `-materi_displacement`, `-materi_stress`, `-condif_temperature`,
`-beam_rotation`, ...). `-none` writes no solution field. Without the record every field is
written. `mesh deform` is an "other" result ([control_print_gid_other](control_print_gid_other.md)),
not a dof field, and stays.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_dof  index  name_0 name_1 ...
control_print_gid_dof  index  -none
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `name_i` | Initialisation name of a field, or `-none`. |

## Example

```
control_print_gid      10  -yes
control_print_gid_dof  10  -materi_displacement -condif_temperature
```

## Differences with the Professional version

None: same selection as the Professional 02-08-2026.

## Tests

`validation-suite/test-2026/`: `control_print_gid_dof1`, `control_print_gid_dof2` (checked by `scripts/build_safe.sh` in a scratch copy).
