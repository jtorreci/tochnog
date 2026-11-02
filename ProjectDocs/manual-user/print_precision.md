# print_precision

## Description

`print_precision number_of_values` (manual Professional 6.997): the number
of significant digits of every printed value (`.dbs`, `control_print`,
GiD, VTK, history files, ...). Example of the manual: with 4 the value
98.123456789 is printed as 98.12. Default 12. The calculation itself is
not affected.

## Usage

```
print_precision 4
```

Test: `validation-suite/test-2026/print_precision1`.

## Differences with the Professional version

The Professional 02-08-2026 prints `number_of_values` decimals in
scientific form (`5.2021e+01`); the GNU prints `number_of_values`
significant digits as in the manual example (`52.02`) (DIVERGENCES
`PRINT-PREC`, class C). The GNU default `.dbs` format (12 significant
digits) is unchanged.
