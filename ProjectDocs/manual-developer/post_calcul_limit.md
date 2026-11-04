# post_calcul_limit

## Implementation (F5 T2 V8)

- **Record:** `database.cc`, DOUBLE, variable length (`DATA_ITEM_SIZE`), no
  index, class POST.
- **Consumer:** `calcul.cc`, `calculate()`: the record is read once; a
  length below `2*nmat` (nmat = dof/operat pairs of `post_calcul`) stops
  with an error and `db_error()`. Per item the file statics
  `calcul_modify_lower / _upper` are set; `calculate_modify()` clamps every
  value of the item (after the multiply, before the absolute value).

## Verification

Professional 02-08-2026 probes: the pair is per dof/operat item (5 values,
3 items: the third pair applies to `size_tot`), the order multiply / limit
/ absolute, and the error for a short record. Tests `post_calcul_modify1`,
`post_calcul_modify2`, `post_calcul_modify_error2`.
