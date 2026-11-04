# post_calcul_absolute

## Implementation (F5 T2 V8)

- **Record:** `database.cc`, INTEGER, length 1, no index, class POST
  (enum `POST_CALCUL_ABSOLUTE` in `tochnog.h` / `tochnog-mod.h`).
- **Consumer:** `calcul.cc`, `calculate()` reads the switch once and sets
  the file static `calcul_modify_absolute` per item. `calculate_modify()`
  applies `fabs()` to every value of the item, after the multiply and the
  limit. It is called right after `calculate_operat()` in the post-point
  loop and in `parallel_calcul_node()`. `refine_l.cc` calls
  `calculate_operat()` directly and gets no modifiers.
- `fabs()` and not `scalar_dabs()`: `scalar_dabs(-0.)` keeps the sign and
  prints `-0`.

## Verification

Professional 02-08-2026 probes (multiply/limit/absolute combinations,
order measured with negative factors and a limit above the value):
`post_calcul_modify1`, `post_calcul_modify2` (`scripts/build_safe.sh`).
