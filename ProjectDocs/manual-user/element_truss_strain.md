# element_truss_strain + element_truss_strain_temperature

## Description

Two print records filled for every truss element after each step:

- `element_truss_strain <element> strain` (manual Professional 6.451):
  the total strain, the length increase divided by the initial length.
- `element_truss_strain_temperature <element> strain` (6.452): the
  accumulated thermal strain, `alpha * dT` with `alpha` from
  `group_truss_expansion` (zero without it or without
  `condif_temperature`).

They are output records; you do not give them in the input.

## Example

`element_truss_strain1.dat` (validation-suite/test-2026): a truss of
length 1, `group_truss_expansion 0 0.001`, the temperature raised by 10
and the free end moved 0.01. Result: `element_truss_strain 1 0.01`,
`element_truss_strain_temperature 1 0.01` (the force is 0: the whole
elongation is thermal).

## Differences with the Professional version

None found (same values on the probe above, Professional 02-08-2026).
