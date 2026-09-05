# print_database_calculation + print_gid_calculation

## Description

Global switches of the output written at the end of a successful
calculation:

- `print_database_calculation` (manual Professional 6.972): `-yes`
  (default) writes the database to `<name>.dbs` after the calculation;
  `-no` skips it.
- `print_gid_calculation`: `-yes` (default) writes the GiD
  (`.flavia`) result files at the end; `-no` skips them.

Useful for large performance runs where the final text dump would
dominate the runtime or the disk usage.

## Usage

```
print_database_calculation [-yes | -no]
print_gid_calculation [-yes | -no]
```

## Example

See `large1.dat`: a 30x30x30 brick performance test that disables both
outputs (the calculation itself remains solver-bound).
