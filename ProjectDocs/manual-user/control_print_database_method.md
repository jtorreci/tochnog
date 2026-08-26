# control_print_database_method

## Description

`control_print_database_method` selects how the database is printed when
a database dump is written (method of the database print). It is useful
to inspect a calculation: `-all` prints every base record, `-size_tot`
reports the memory size of every base record (plus the size of the
system matrix), and `-size_tot_large` reports only the records larger
than 1 Mb (plus the size of the system matrix).

If the record is not specified, the method defaults to `-all`.

## Uso

Place it in the data part with the same `icontrol` index as the
`control_timestep` record:

```
control_print_database_method 20 -size_tot
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. |
| method    | `-all` (default): print all database base records. `-size_tot`: print the size of all database base records. `-size_tot_large`: print the size of the database base records larger than 1 Mb. |

With `-size_tot` and `-size_tot_large` the size of the system matrix
(number of equations of the last solve) is also printed.

## Output

For `-all` the file `<inputfile><index>.dbs` contains all base records
(the complete status of the calculation). For `-size_tot` /
`-size_tot_large` the same file contains one `Size of <record> is
<bytes>` line per record (only records > 1 Mb for `-size_tot_large`),
a `Size of the system matrix is <n>` line and a `Total size is <n>`
footer.

## Example

```
control_print_database_method 20 -size_tot
control_timestep                20  0.1 0.1
```

Produces `dbmeth20.dbs` with the memory footprint of every base record
of the model.

## Validation

Test `dbmeth` (validation-suite/test-2014): the three methods are
exercised on the same small model and the generated `.dbs` files are
checked: `-all` prints every record without size lines; `-size_tot`
prints the sizes plus the system matrix; `-size_tot_large` prints only
the system matrix line (all records of the model are smaller than
1 Mb — the A/B against `-size_tot` discriminates the size filter).
