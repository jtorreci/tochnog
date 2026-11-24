# area_element_group family (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — new enums
  `AREA_ELEMENT_GROUP_ELEMENT/_INTERFACE/_NODE/_TIME`,
  `AREA_ELEMENT_GROUP_SEQUENCE_ELEMENT_GROUP` (Professional alias),
  `AREA_ELEMENT_GROUP_SEQUENCE_GEOMETRY_METHOD`,
  `AREA_ELEMENT_GROUP_SEQUENCE_INTERFACE`; switch enum
  `ANY_BUT_NOT_ALL`.
- `database.cc` — registrations. NOTE: `AREA_ELEMENT_GROUP_METHOD` and
  `AREA_ELEMENT_GROUP_SEQUENCE_METHOD` were READ by group.cc since the
  GNU origins but NEVER registered — nobody could use them (fixed
  here). The `data_required` of SEQUENCE_GEOMETRY/TIME points to the
  legacy ELEMENTGROUP; since the Professional alias is a separate
  record the hard check would reject alias-only inputs, so the
  combination is validated at runtime (group.cc has its own "should be
  specified" error).
- `group.cc`:
  - `area_element_group()`: `_element` name filter (el[0]), `_node`
    list (array_member instead of geometry()), `_method` extended with
    `-any_but_not_all` and positive-integer-N (count_in), `_interface`
    guard (skip elements whose current group has GROUP_INTERFACE when
    the switch is not -yes).
  - new `area_element_group_time_active()`: any `_time -yes` record?
  - `area_element_group_sequence()`: the Professional alias copy runs
    BEFORE the legacy max_index query (an alias-only input would
    otherwise see max=-1 and never enter the loop — found by
    debugging); `_geometry_method` overrides the method for the
    geometry test; methods extended as above; `_interface` guard.
- `top.cc` — `step_close`: when `area_element_group_time_active()`,
  re-run `area_element_group(VERSION_NORMAL)` per step.
- `bounda.cc` — `bounda_time_factor` (manual Professional 6.36):
  multiplies the LOAD values (odd positions; times untouched) of the
  `bounda_time` record with the same index, after `_units`. Handles the
  single-value constant form too.

## Gotchas found (test-side, documented for future tests)

- `bounda_unknown -ra a b` is a RANGE (a..b inclusive), not a pair of
  nodes. For non-contiguous nodes use several `-ra a a` ranges or
  several records. A stray range silently prescribes dofs on unintended
  nodes (strain=0 mysteries).
- With the Z node convention (see
  condif_convection_radiation_edge_normal.md) the "top" edge of quad4
  1,2,3,4 is nodes 3,4; the base is 1,2.

## Verification (suite 79/79)

- `aeg_node`: `_node`-less variant with `_method -any` + `_time -yes`:
  brick covers the left column only -> element switches to the group
  with young 2000: sigyy -4.0 vs -2.0 (the A/B proves both the
  regrouping and the per-step re-evaluation).
- `aeg_seq`: Professional alias `..._element_group` + `_geometry_method
  -any`; switch at t=0.1 with 3 steps: sigyy = -(2 steps E1000 + 1 step
  E2000) = -2.0 exact (also proves the alias copy order fix).
- `bt_factor`: table load 5.0 with factor 2.0 -> effective 10 per node:
  sigyy 20.0 (vs 10 without the factor).
