# control_timestep_iterations_automatic (3 values)

## Implementation

- **Registration** (database.cc): `CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC`
  data_length changed from 2 to 3 (manual Professional 6.386:
  ratio_criterium minimal_timestep maximum_timestep).
- **Consumption** (top.cc, CONTROL_TIMESTEP branch): the record is read
  into `control_timestep_iterations_automatic[]` and the slots are
  mapped as
  ```c
  ratio_criterium = control_timestep_iterations_automatic[0];
  maximum_timestep = control_timestep_iterations_automatic[2];
  ```
  (`minimal_timestep` at [1] is not consumed by the current automatic
  scheme, which only reduces the step; the reduction uses the ratio
  criterion and the maximum bound).

## Verification

- slope_classical_numerical: with the corrected slot mapping the model
  advances to completion (previously the maximum was read as 1.e-6 and
  the run timed out at 45 s). The collapse time target (1.275) is not
  met (measured 2.0) - the phi-c reduction physics of the slope model
  still needs calibration (documented in SEGUIMIENTO).
- dam_building: the layered construction uses
  `control_timestep_iterations_automatic` with 3 values on every layer
  (parsed and run).

## Pending

- `minimal_timestep` (slot [1]) is registered but not consumed: the GNU
  automatic scheme reduces the step without a lower bound. Documented
  partial.
