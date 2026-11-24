# control_groundflow_consolidation_apply

## Files and functions

- `groundfl.cc` — in `groundflow()` (lines 62-68). Read with the timestep
  control index:
  ```c
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_GROUNDFLOW_CONSOLIDATION_APPLY, icontrol,
    &control_groundflow_consolidation_apply, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  ...
  if ( control_groundflow_consolidation_apply==-NO ) materidivergence = -NO;
  ```
- `database.cc` — keyword registration (alphabetical, between
  `CONTROL_DISTRIBUTE_VALUES` and `CONTROL_EIGEN`): type `INTEGER`,
  `data_length = 1`, `data_class = CONTROL`,
  `data_required = GROUNDFLOW`.
- Enum `CONTROL_GROUNDFLOW_CONSOLIDATION_APPLY` in `tochnog.h` (alphabetical,
  same position as the registration).

## Implementation details

- Indexed by the timestep control index (`icontrol`), so different timesteps
  can switch consolidation independently.
- `GET_IF_EXISTS`: absent record keeps the default `-NO`, so the global
  `groundflow_consolidation_apply` (or the default `-NO`) governs. The default
  is `-NO` since the u-p consolidation sprint (commit `0560dcf`); the GNU
  legacy initializer was `-YES`.
- Same precedence chain as `groundflow_consolidation_apply`; this record sits
  above the global one but below the legacy `CONTROL_OPTIONS_SKIP_...`.

## External dependencies

- Core `db()` accessor and the `ICONTROL` global (set by the timestep loop).

## Hardcoded parameters / pending refactorings

- Default `-NO` is implicit in the variable initializer in `groundflow()`.
