# post_calcul_multiply

## Implementation (F5 T2 V8)

- **Record:** `database.cc`, DOUBLE, variable length (`DATA_ITEM_SIZE`), no
  index, class POST.
- **Consumer:** `calcul.cc`, `calculate()`: read once, a length below
  `nmat` stops with an error and `db_error()`. Per item the file static
  `calcul_modify_factor` is set; `calculate_modify()` multiplies every
  value of the item first.

## Verification

Professional 02-08-2026 probes: 3 and 5 factors for 3 items give the same
result (per dof/operat item, extra values ignored); 2 factors are
rejected; a negative factor with `post_calcul_absolute` shows the absolute
value is taken after the factor. Tests `post_calcul_modify1`,
`post_calcul_modify_error1`.
