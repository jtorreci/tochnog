# control_timestep_until_data — developer

- Records `CONTROL_TIMESTEP_UNTIL_DATA` (INTEGER, up to 30 values),
  `CONTROL_TIMESTEP_UNTIL_MINIMUM` / `_MAXIMUM` (DOUBLE, up to 10), all
  CONTROL and `fixed_length = 0` (`database.cc`).
- `control_timestep_until_reached( icontrol )` (`top.cc`), called right
  after `step_close( YES, ... )` of a converged step; 1 breaks the
  timestep `while` loop of the index. `db_error` when a minimum / maximum
  record is missing or its length differs from the number of triples.
- Values through `data_item_value()` (`bounda.cc`), factored out of
  `bounda_time_until_data` (same lookup: dof label through `dof_label`,
  principal slot for a record of `npuknwn` values, or a plain position).
