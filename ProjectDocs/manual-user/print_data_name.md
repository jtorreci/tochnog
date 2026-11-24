# print_data_name

## Description

`print_data_name switch` (manual Professional 6.971): with `-yes` every
data name of the database is written to `tochnog_data_name.txt` (an empty
line, `Names:`, then one name per line, in the database order). Search it
for options, e.g. `grep group tochnog_data_name.txt`. The list also holds
internal names that are not input records; check this manual.

## Usage

```
print_data_name -yes
```

Test: `validation-suite/test-2026/print_data_name1`.

## Differences with the Professional version

Same layout; the names are the GNU database names (about 1650 in F5, the
Professional 02-08-2026 lists 1591).
