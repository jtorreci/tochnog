# group_groundflow_nonsaturated_eps_permeability

## Files and functions

- `groundda.cc` — in `groundflow_data()` (lines 101-104). Read when present
  and applied as a floor to the van Genuchten relative permeability:
  ```c
  if ( db_active_index( GROUP_GROUNDFLOW_NONSATURATED_EPS_PERMEABILITY,
      gr, VERSION_NORMAL ) )
    get_group_data( GROUP_GROUNDFLOW_NONSATURATED_EPS_PERMEABILITY, gr, element,
      new_unknowns, &eps_permeability, ldum, GET_AND_CHECK );
  ...
  if ( krel < eps_permeability ) krel = eps_permeability;
  ```
- `database.cc` — keyword registration: type `DOUBLE_PRECISION`,
  `data_length = 1`, `data_class = GROUNDFLOW`,
  `data_required = GROUP_TYPE`.
- `check.cc` — requires `groundflow_pressure`.
- Enum `GROUP_GROUNDFLOW_NONSATURATED_EPS_PERMEABILITY` in `tochnog.h`.

## Implementation details

- Only meaningful together with `group_groundflow_nonsaturated_vangenuchten`;
  without it the value is read but unused.
- Floors `krel` AFTER the Mualem expression and AFTER the `[0,...]` clamp, so
  it is the effective lower bound of the relative permeability.

## External dependencies

- Core `db()` accessor, `get_group_data()`.

## Hardcoded parameters / pending refactorings

- Default (record absent) is `eps_permeability = 0` (no floor), so an absent
  record behaves as an unbounded reduction.
