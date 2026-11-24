# post_calcul_safety_default

## Implementation (F5 T2 V8)

- **Record:** `database.cc`, DOUBLE, fixed length 2 (`fixed_length` 1), no
  index, class POST.
- **Consumer:** `calcul.cc`, `calculate_operat()`, safety branch: the
  defaults `TINY` / 0 are replaced by the record; the test is
  `|p_div| < eps` (p_dynamic for piping, p_total for lifting).

## Verification

`post_calcul_safety_default1`: node 1 and node 41 `node_dof_calcul` equal
the Professional 02-08-2026 (`scripts/build_safe.sh`).
