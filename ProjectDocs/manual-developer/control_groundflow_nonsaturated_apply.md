# control_groundflow_nonsaturated_apply

## Files and functions

- `groundda.cc` — in `groundflow_data()` (lines 91-93). Read with the timestep
  control index as the per-timestep gate of the van Genuchten law:
  ```c
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_GROUNDFLOW_NONSATURATED_APPLY, icontrol, &nonsaturated_apply,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  ```
- `database.cc` — keyword registration (alphabetical, between
  `CONTROL_DISTRIBUTE_VALUES` and `CONTROL_EIGEN`): type `INTEGER`,
  `data_length = 1`, `data_class = CONTROL`,
  `data_required = GROUNDFLOW`.
- Enum `CONTROL_GROUNDFLOW_NONSATURATED_APPLY` in `tochnog.h`.

## Implementation details

- Indexed by the timestep control index (`icontrol`), so different timesteps
  can switch the non-saturated model independently.
- `GET_IF_EXISTS`: absent record keeps the default `-YES`, so the global
  `groundflow_nonsaturated_apply` (or the default) governs.
- Applied as an OR-gate together with the global switch: the law runs only if
  BOTH end up `-YES`.

## External dependencies

- Core `db()` accessor and the `ICONTROL` global.

## Hardcoded parameters / pending refactorings

- Default `-YES` is implicit in the variable initializer in
  `groundflow_data()`.
