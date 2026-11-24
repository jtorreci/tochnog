# element_truss_strain + element_truss_strain_temperature (developer)

## Where

- `truss.cc`, `truss()`: after `ELEMENT_TRUSS_FORCE` is PUT (VERSION_NEW),
  `ELEMENT_TRUSS_STRAIN = (|x2-x1|_new - L0) / L0` and
  `ELEMENT_TRUSS_STRAIN_TEMPERATURE = old value (VERSION_NORMAL) +
  group_truss_expansion * inc_temperature` (the step temperature increment
  already computed for the thermal stretch).
- `top.cc`: `db_allocate( ..., VERSION_NEW, MINIMAL )` next to
  `ELEMENT_TRUSS_FORCE` (no allocation inside the parallel element loop).
- `database.cc`: DOUBLE_PRECISION, length 1, `version_all`, class ELEMENT;
  the step-end version copy promotes them like the truss force.
- Enums appended before `LAST_DUMMY` in `tochnog.h` / `tochnog-mod.h`.

## Notes

- L0 is the initial length of the truss (`initial_coord`), the same for all
  memory types.
- Verified against the Professional 02-08-2026 (`scratchpad/f5/v4` probe
  `truss`); test `element_truss_strain1`.
