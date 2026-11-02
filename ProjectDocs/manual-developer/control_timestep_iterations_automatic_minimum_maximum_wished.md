# control_timestep_iterations_automatic_minimum_maximum_wished — developer

- Record `CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC_MINIMUM_MAXIMUM_WISHED`
  (INTEGER, 3 values, fixed length, CONTROL; `database.cc` after
  `CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC_STOP`).
- `top()` (`top.cc`), automatic branch: `iteration_min`, `iteration_max`,
  `iteration_wished` default 2 8 4 (was 2 20 4 and the `ITERATION_WISHED`
  macro, removed); the record overrides them; `db_error` when maximum < 2,
  minimum < 1 or minimum > maximum. The step factor uses
  `iteration_wished`.
- `CONTROL_TIMESTEP_ADJUST_MINIMUM_ITERATIONS` (INTEGER, CONTROL) is
  registered and not consumed.

## Pending (BACKLOG C94)

- The Professional step-size law (x1.25 / x1.5625 growth, /4 after the
  maximum) is not reproduced; `control_timestep_adjust_minimum_iterations
  -yes` (Professional default) has no specified rule.
