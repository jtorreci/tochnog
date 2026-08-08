# bounda_dof

## Files and functions

- `bounda.cc` — within `bounda()`: `bounda_dof` is implemented as a functional
  alias of `bounda_unknown`.
  - Local `max_bounda_dof=0` (line 29), added to the `max_bounda` calculation:
    `db_max_index( BOUNDA_DOF, max_bounda_dof, ... )` (line 64) and
    `if ( max_bounda_dof>max_bounda ) max_bounda = max_bounda_dof;` (line 83).
  - The iboun loop treats both records as "unknown" (Dirichlet):
    `unknown = db_active_index( BOUNDA_UNKNOWN, iboun, ... ) ||
               db_active_index( BOUNDA_DOF, iboun, ... );` (lines 86–87).
  - Value read (lines 183–188): uses `BOUNDA_DOF` when it is active, otherwise
    falls back to `BOUNDA_UNKNOWN`.
- `database.cc` — keyword registration (lines 188–192): name `bounda_dof`,
  `type = INTEGER`, `data_length = MBOUNDA`, `fixed_length = 0`,
  `data_class = BOUNDA`.
- `tochnog.h` — enum `BOUNDA_DOF` (line 142), before `BOUNDA_UNKNOWN` (line 154).
- `tochnog-mod.h` — mirror enum `BOUNDA_DOF` (line 135), must stay in sync.

## Implementation details

- `bounda_dof` is a pure alias: the record is read into the same `val` array
  and processed by the identical code path as `bounda_unknown`, so it supports
  the same selection modes — node ranges (`-range`), `-all`, `-node_set` and
  geometries (`-geometry_line`, `-geometry_point`, ...) — through the shared
  `val[0]` dispatch in the node loop (lines 288–305).
- The `unknown` flag drives both the value read and the dof-label lookup; only
  one of `bounda_dof`/`bounda_unknown` or `bounda_force` may be active per
  boundary index.
- The dofs prescribed are the primary dofs (e.g. `-velx`, `-temp`, `-pres`)
  resolved through `dof_label`/`array_member`, identical to `bounda_unknown`.
- Values over time come from the matching `bounda_time` record, as for
  `bounda_unknown`.

## External dependencies

None. Reuses the existing `bounda_unknown` code path and the `geometry()`
routine for geometry-based selection.

## Hardcoded parameters / pending refactorings

- PENDING: error reporting still references the GNU name — the check
  `if ( bounda_length<2 ) db_error( BOUNDA_UNKNOWN, iboun );` (line 189) and
  the dof lookup error `db_error( BOUNDA_UNKNOWN, iboun );` (line 355) fire
  even when the input used `bounda_dof`. Messages may confuse users of the
  professional keyword.
- The alias duplicates the `bounda_unknown` registration pattern in
  `database.cc`; the two records could share a single registration helper to
  avoid drift between `BOUNDA_DOF` and `BOUNDA_UNKNOWN`.
- `data_required` is unset for `BOUNDA_DOF`, so `bounda_dof` works standalone
  (like `bounda_unknown`); keep `BOUNDA_CONSTANT`'s
  `data_required = BOUNDA_UNKNOWN` (database.cc:204) in mind if records ever
  need to require either name.
