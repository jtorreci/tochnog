# groundflow_nonsaturated_apply

## Files and functions

- `groundda.cc` — in `groundflow_data()` (lines 88-93). Global gate of the
  van Genuchten law:
  ```c
  db( GROUNDFLOW_NONSATURATED_APPLY, 0, &nonsaturated_apply, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_GROUNDFLOW_NONSATURATED_APPLY, icontrol, &nonsaturated_apply,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  if ( nonsaturated_apply==-YES ) { ... law ... }
  ```
- `database.cc` — keyword registration: type `INTEGER`, `data_length = 1`,
  `no_index = 1`, `data_class = GROUNDFLOW`.
- `check.cc` — requires `groundflow_pressure`.
- Enum `GROUNDFLOW_NONSATURATED_APPLY` in `tochnog.h`.

## Implementation details

- `no_index`: read with the fixed index `0`.
- `GET_IF_EXISTS`: absent record keeps the default `-YES`.
- The per-timestep `CONTROL_GROUNDFLOW_NONSATURATED_APPLY` overrides the global
  value when present; otherwise the global (or the default) applies.

## External dependencies

- Core `db()` accessor; the `ICONTROL` global (set by the timestep loop).

## Hardcoded parameters / pending refactorings

- Default `-YES` is implicit in the variable initializer in
  `groundflow_data()`.
