# post_calcul_materi_stress_force_element_group

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`, one
  record shared by all `-force` post_calcul records): INTEGER,
  `data_length = DATA_ITEM_SIZE`, `fixed_length = 0` (a list of element
  groups), `data_required = POST_CALCUL`.
- Mandatory record: `post_calcul_materi_stress_force_validate()`
  (`calcul_force.cc`, called once per `-force` record in `calculate()`)
  aborts with a clear message when the record is absent or lists no
  group.
- Consumed by the numerical integration in lots 2/3 (element sides of
  the listed groups); lot 1 only validates and stores.

## Pendiente

- The actual element-side selection over the listed groups lands in
  lot 2/3.
