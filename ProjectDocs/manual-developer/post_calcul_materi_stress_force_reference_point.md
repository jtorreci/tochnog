# post_calcul_materi_stress_force_reference_point

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`):
  DOUBLE_PRECISION, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0` (variable: `ngroups * ndim` values, one point per
  element group).
- Validated in `post_calcul_materi_stress_force_validate()`
  (`calcul_force.cc`): total length must be `ngroups * ndim` (clear
  error otherwise); in 3D the record is REQUIRED (error when absent),
  in 2D a warning is issued and the documented default (0,0) is used.
- Consumed in lots 2/3 to orient the forces/moments consistently
  outwards/inwards in the thickness direction (the far/outer nodes of
  [`post_calcul_materi_stress_force_outer`](post_calcul_materi_stress_force_outer.md)
  are also selected with it).

## Pendiente

- The orientation logic (sign of the vectors relative to the reference
  point) lands in lot 2/3.
