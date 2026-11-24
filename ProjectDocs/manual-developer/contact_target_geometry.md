# contact_target_geometry / contact_target_geometry_switch (developer)

## Where

- `contact.cc`, `static long int contact_geometry_get(itar, geom[], swit)`:
  returns the target geometry of index `itar` and its switch, from
  `CONTACT_GEOMETRY`/`CONTACT_GEOMETRY_SWITCH` (GNU name, wins) or from
  `CONTACT_TARGET_GEOMETRY`/`CONTACT_TARGET_GEOMETRY_SWITCH` (Professional
  name, manual 6.106/6.107). Read-only (`GET` / `GET_IF_EXISTS`).
- `contact.cc`, `parallel_contact()`: the geometry loop runs up to the
  larger of the two `db_max_index` values and calls `contact_geometry_get`
  for the target test; the switch read inside the `inside_target` branch
  uses the value it returned.
- Records: `database.cc` (`CONTACT_TARGET_GEOMETRY`, `INTEGER`, length 2;
  `CONTACT_TARGET_GEOMETRY_SWITCH`, length 1, `data_required` on the
  target geometry).

## Implementation details

The old alias (added with the `contact_*` family, `23262bf`) copied
`CONTACT_TARGET_GEOMETRY 0` to `CONTACT_GEOMETRY 0` with a db `PUT` at the
start of `parallel_contact()`. That routine runs through
`parallel_sys_routine` (top.cc), where `parallel_active` is set, so the
first `PUT` allocated and `db_allocate()` aborted with "Data is allocated
in a parallel loop" (6 audit runs: `vs/contact*`, `examp6/9/24` with the
renamed record). Indices > 0 were never copied (validation_10 and
pressed_rubber_cylinder lost their contact line; contac2 its point).

The new lookup reads both names at the point of use, so no record is
created and the dbs output keeps the user's spelling. Same pattern as
`contact_heat_generation` (dual read, no alias PUT).

## Verification

- test-2026 `contact_target_geometry_index1` (index 1, x = 1.0),
  `contact_target_geometry_index0` (index 0, x = 1.0),
  `contact_target_geometry_switch1` (index 0 + switch, x = 1.02): RED on
  the old binary (rc 1: free flight / parallel-loop abort), GREEN on the
  new one, Professional 02-08-2026 meets all three targets.
- Rename neutrality: `contact_geometry` -> `contact_target_geometry` in
  `contact`, `contact_block`, `contact_ctrl_apply`, `contact_heatgen`,
  `examp6`, `examp9` gives byte-identical `node_dof` records.
- Pro inputs: `validation_10` sigyy -21232 (Pro -21031, target -20538 +-
  1e4; before ~0), `pressed_rubber_cylinder` max |dof - Pro| 6e-7 (before
  0.10), `contac2` x = 1.0 (= Pro).

## Pending

- Element targets with `contact_target_element_group` need the mesh
  boundary: see `contact_target_element_group.md`.
