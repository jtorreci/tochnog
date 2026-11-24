# group_materi_elasti_camclay_pressure_min

## Implementación

- **Record**: registered in `database.cc`
  (`GROUP_MATERI_ELASTI_CAMCLAY_PRESSURE_MIN`: DOUBLE_PRECISION,
  length 1, data_class MATERI, data_required GROUP_TYPE), right after
  `group_materi_elasti_camclay_poisson`.
- **Hook**: the camclay elastic block in `set_stress()` (stress.cc,
  the block gated by `GROUP_MATERI_ELASTI_CAMCLAY_G` /
  `GROUP_MATERI_ELASTI_CAMCLAY_POISSON`). Before computing
  `k = (1.+e)*pressure/kappa`, if the record exists the pressure is
  clamped: `if ( pressure<pressure_min ) pressure = pressure_min;`.
  The clamp applies to the pressure used for K only; the stress state
  itself is unchanged.
- **check.cc**: `GROUP_MATERI_ELASTI_CAMCLAY_PRESSURE_MIN` requires
  `materi_stress` and `materi_history_variables` (same state as the
  camclay elastic records).

## Física

Manual (6.647): "This specifies a minimal allowed value for the
pressure in the calculation of the bulk modulus for the camclay model.
In the calculation pressures below pressure_min will be set to
pressure_min. This prevents numerical problems for very low bulk
modulus K values." The camclay bulk modulus K = (1+e)*p/kappa is
pressure-proportional: as p -> 0, K -> 0 (near-singular C, poisson ->
-1); for tensile states K would turn negative. The clamp keeps K >=
(1+e)*pressure_min/kappa.

## Gotchas

- The camclay elastic block derives poisson from (k, g) and builds C
  from (young, poisson); the resulting bulk modulus of C equals k
  exactly (consistent conversion). The clamp keeps that k finite.
- Without the record the degenerate response is stable but physically
  wrong (poisson ~ -1, sign-flipped lateral stress); the A/B pair
  mc_pressure_min / mc_pressure_min_off discriminates.
