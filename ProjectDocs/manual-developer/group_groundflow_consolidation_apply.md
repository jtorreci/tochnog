# group_groundflow_consolidation_apply

## Files and functions

- `groundfl.cc` — in `groundflow()` (lines 51-60). Read per element group next
  to the existing `GROUP_GROUNDFLOW_MATERIDIVERGENCE`:
  ```c
  db( GROUP_GROUNDFLOW_CONSOLIDATION_APPLY, gr, &group_groundflow_consolidation_apply,
    ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  ...
  if ( group_groundflow_consolidation_apply==-NO ) materidivergence = -NO;
  ```
  `gr` is the element group of the current integration point.
- `database.cc` — keyword registration (alphabetical, between
  `GROUP_GROUNDFLOW_CAPACITY_NONLINEAR_PARAMETERS` and
  `GROUP_GROUNDFLOW_MATERIDIVERGENCE`): type `INTEGER`, `data_length = 1`,
  `data_class = GROUNDFLOW`, `data_required = GROUP_TYPE`.
- `check.cc` — requires `groundflow_pressure` and `materi_velocity`.
- Enum `GROUP_GROUNDFLOW_CONSOLIDATION_APPLY` in `tochnog.h`.

## Implementation details

- Indexed by element group (`gr`); the group resolution follows the standard
  `group_type ... -groundflow` mechanism.
- `GET_IF_EXISTS`: absent record keeps default `-NO` (since the u-p
  consolidation sprint, commit `0560dcf`; the GNU legacy initializer was
  `-YES`, aligned with the Professional default `-no`).
- Lowest precedence of the three consolidation switches: if the global or
  timestep switch sets `-NO`, the group cannot re-enable it (only `-NO` is
  propagated to `materidivergence`).

## External dependencies

- Core `db()` accessor only; `gr` comes from the element assembly loop.

## Hardcoded parameters / pending refactorings

- Mirrors `GROUP_GROUNDFLOW_MATERIDIVERGENCE` with inverted semantics; the two
  could be unified in a single switch to avoid confusion.
