# post_calcul_label

## Description

`post_calcul_label label_0 label_1 ...` (manual Professional 6.904). An
output record: after a calculation the database file (`.dbs`) holds the
names of the items calculated by `post_calcul`, one name per value, in the
order of `node_dof_calcul` / `post_point_dof_calcul`.

You do not give this record yourself. When an input file holds it (for
example a database file read back), its names are ignored.

## Usage

The record appears in the `.dbs` file, for example:

```
post_calcul_label  -aept -va0ept -va1ept -va2ept -sept
```

for `post_calcul -materi_strain_total -average -materi_strain_total -prival
-materi_strain_total -size_tot` in 1D. The names can be used as
`target_item` labels and in `control_print_vtk_dof_calcul`.

## Example

`validation-suite/test-2026/post_calcul_modify1.dat` (the input holds a
`post_calcul_label` copied from the Professional, which is ignored; the
`.dbs` gets the GNU names).

## Differences with the Professional version

The generic operators are named without an underscore, as in the manual
example (`-aept`). The Professional 02-08-2026 writes `-a_ept -va0_ept
-s_ept`. The pressure, safety and apparent names are the same (`-to_pres`,
`-safety_piping`, `-young_apparent`). DIVERGENCES POST-CALCUL-LABEL-NAMES, BACKLOG C109.
