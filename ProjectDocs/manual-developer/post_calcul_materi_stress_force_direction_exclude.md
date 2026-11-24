# post_calcul_materi_stress_force_direction_exclude

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`):
  DOUBLE_PRECISION, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0` (variable: exactly `ndim` values in 3D).
- Validated in `post_calcul_materi_stress_force_validate()`
  (`calcul_force.cc`): length must be `ndim`; in 2D the direction
  records warn and are ignored; exclude XOR include with
  [`post_calcul_materi_stress_force_direction_include`](post_calcul_materi_stress_force_direction_include.md)
  (manual 6.913 "not both").
- Semantics (consumed in lots 2/3): exclude the element sides with
  `|n . dir| > 1 - eps`, `eps` default `1.e-8`
  (`post_calcul_materi_stress_force_direction_exclude_epsilon`).

## Pendiente

- Side normal computation and the `|n . dir|` test land in lot 2/3.
