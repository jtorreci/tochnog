# groundflow_consolidation_apply

## Files and functions

- `groundfl.cc` — in `groundflow()` (lines 49-67). The switch is read after the
  existing `GROUP_GROUNDFLOW_MATERIDIVERGENCE` and applied on top of it:
  ```c
  db( GROUNDFLOW_CONSOLIDATION_APPLY, 0, &groundflow_consolidation_apply,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  ...
  if ( groundflow_consolidation_apply==-NO ) materidivergence = -NO;
  ```
  `materidivergence==-NO` skips the divergence block (lines 83-100), which is
  the consolidation coupling between `materi_velocity` divergence and the
  pressure equation.
- `database.cc` — keyword registration (after `GROUNDFLOW_ADDTOPRESSURE`):
  type `INTEGER`, `data_length = 1`, `no_index = 1`,
  `data_class = GROUNDFLOW`.
- `check.cc` — requires `groundflow_pressure` to be initialized.
- Enum `GROUNDFLOW_CONSOLIDATION_APPLY` in `tochnog.h` (alphabetical, between
  `GROUNDFLOW_ADDTOPRESSURE` and `GROUNDFLOW_DENSITY`).

## Implementation details

- `no_index`: the database accessor is called with the fixed index `0`.
- Read with `GET_IF_EXISTS`: when the record is absent, the variable keeps its
  default `-YES` and nothing changes.
- Precedence (from lowest to highest): `GROUP_GROUNDFLOW_MATERIDIVERGENCE`,
  `GROUP_GROUNDFLOW_CONSOLIDATION_APPLY` (group), `GROUNDFLOW_CONSOLIDATION_APPLY`
  (global), `CONTROL_GROUNDFLOW_CONSOLIDATION_APPLY` (timestep),
  `OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE` /
  `CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE` (legacy skip switches).
  Any `-NO`/`-YES`-skip forces `materidivergence = -NO`.

## External dependencies

- Core `db()` accessor only; no external library.
- Globals `materidivergence`, `materi_velocity`, `pres_indx` (in `tochnog.h`).

## Hardcoded parameters / pending refactorings

- The three switches (group/global/timestep) plus the two legacy skip switches
  all funnel into the single `materidivergence` flag; the precedence chain is
  simple but implicit. A small helper would make the intent clearer.
