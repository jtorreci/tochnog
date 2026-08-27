# post_calcul_materi_stress_force_direction_exclude_epsilon

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`):
  DOUBLE_PRECISION, `data_length = 1` (fixed, exactly one value).
- Default `1.e-8` (manual 6.910) applies when the record is absent;
  the default lives in the documentation of the exclusion test
  (lot 2/3), the record is validated and stored in lot 1.
