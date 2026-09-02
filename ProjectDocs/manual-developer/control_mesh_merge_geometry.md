# control_mesh_merge_geometry family

## Where implemented

- `database.cc` (db_initialize): `control_mesh_merge_geometry_not`
  (INTEGER, variable length, pairs of geometry entity + index) is the
  Professional name (manual 6.215) of the record previously registered
  as `control_mesh_merge_not` (legacy spelling kept as a `db_number`
  alias); new enum `CONTROL_MESH_MERGE_GEOMETRY` for
  `control_mesh_merge_geometry` (manual 6.214, INTEGER, variable). The
  Professional `control_mesh_merge_eps_coord` is aliased to the
  existing `CONTROL_MESH_MERGE_EPSCOORD` in `db_number`.
- `merge.cc` (`merge`): both geometry scopes are evaluated per node
  with `geometry()` (NODE_START_REFINED, PROJECT_EXACT) before the
  coordinate merging; `node_merge_not[]` is set for the blacklist
  (`_geometry_not`) and for nodes OFF the whitelist
  (`_geometry`); the existing `_macro_generate` restriction is
  unchanged and applies afterwards.

## Implementation details

- Both records support multiple entity pairs in one record (loop over
  pairs, not a fixed geometry_entity[2] as before).
- Order of evaluation: blacklist first, whitelist second (only nodes
  still mergeable), macro_generate restriction last.
- `mpc_element_group` and `control_mesh_merge_geometry_not` combine for
  interface problems: nodes on the excluded geometry keep double
  numbers and are then tied with mpc records.

## Verification

- `ground19_water_under_dam`: the `control_mesh_merge_geometry_not 60
  -wall` record parses and the wall nodes stay unmerged (the flow
  field around the wall reproduces the Professional flux; the mpc ties
  make the mechanical field follow).
- Legacy merge inputs (whitelist `control_mesh_merge_geometry` syntax
  of the sfnet suite) parse through the new records; merge2/tutorial_3
  keep their pre-existing blockers (reported in SEGUIMIENTO).

## Pending

- None for the records above; `control_mesh_merge_eps_coord` is
  parse-level only (the underlying eps machinery reads
  `CONTROL_MESH_MERGE_EPSCOORD` per control index; no corpus test
  exercises the Professional spelling).
