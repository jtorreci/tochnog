# post_calcul_materi_stress_force_plot_switch

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`):
  INTEGER, `data_length = DATA_ITEM_SIZE`, `fixed_length = 0`
  (variable: one switch per VECTOR item).
- Validated in `post_calcul_materi_stress_force_validate()`
  (`calcul_force.cc`): length must be 3 (2D) or 4 (3D) - the number of
  VECTOR items (nor, she, mom / nor, she, mom1, mom2), NOT the number
  of result items (9/16); every value must be `-yes`/`-no`.
- Semantics (consumed in lots 2/3): `-yes` inverts the plot vector
  direction of the corresponding item (manual 6.916).

## Implementación (lote 2, 2D)

- Implemented in `msf_element_contribution_2d` (calcul_force.cc):
  `-yes` inverts the x/y plot components of the item (the drawing
  direction, manual 6.916); the s (size) component is untouched.
  Pending for 3D (lot 3).
