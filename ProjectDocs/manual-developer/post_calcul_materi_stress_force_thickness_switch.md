# post_calcul_materi_stress_force_thickness_switch

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`):
  INTEGER, `data_length = DATA_ITEM_SIZE`, `fixed_length = 0`
  (variable: one switch per element group).
- Validated in `post_calcul_materi_stress_force_validate()`
  (`calcul_force.cc`): length must equal `ngroups`; every value must be
  `-yes` or `-no` (`db_error` otherwise).
- Semantics (consumed in lots 2/3): `-yes` selects the LONGEST element
  direction of the evaluated side as structure thickness (manual
  6.917); default `-no` (shortest).

## Pendiente

- The shortest/longest direction analysis lands in lot 2/3.
