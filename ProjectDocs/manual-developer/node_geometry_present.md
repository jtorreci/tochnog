# node_geometry_present + print_node_geometry_present

## Implementation

- Keywords registered in `database.cc`:
  - `PRINT_NODE_GEOMETRY_PRESENT` (INTEGER, no_index, global switch).
  - `PRINT_NODE_GEOMETRY_PRESENT_NODE_TYPE` (INTEGER, no_index, global
    default node_type).
  - `NODE_GEOMETRY_PRESENT` (INTEGER, variable length per node,
    version_all=1, class NODE) — the stored (name, index) pairs.
  Enums appended in `tochnog.h`/`tochnog-mod.h` (same enum, same
  order).
- Fill routine `node_geometry_present_calculate()` (`geometry.cc`),
  hooked at the TOP of `step_start()` (`top.cc`): runs once per step
  when the switch is `-yes`, at a moment where the node coordinates
  still hold the converged state of the previous step. For every node
  and every active geometry entity (the 14 concrete geometry types) it
  calls `geometry()` with the default node_type of the print record
  (`-node_start_refined`); the per-geometry `geometry_node_type`/
  `geometry_projection_type` overrides apply inside `geometry()`. The
  resulting (entity value, index) pairs overwrite the record of the
  node; an empty result deletes it (both behaviors measured on the
  Professional).
- Target access: `-node_geometry_present <node> <position>` reads the
  integer at that position of the node record (position 1 = index of
  the first present geometry). No change to the target machinery was
  needed.
- Semantics calibrated against the Professional binary 25-10-2023 with
  a time sweep (t_end = 3/5/8/9/10): the fill uses the coordinates of
  the start of the last step, and the last fill wins — reproduced by
  filling at step_start (the GNU moves the NODE coordinates at the end
  of each step, so at step_start k they hold the positions of t=k-1).
  A static model (no motion) leaves the record empty exactly like the
  Professional (default `-node_start_refined` coordinates are only
  meaningful after refinement/motion bookkeeping — measured).
- Verified: `node_type_1.dat` rc=0 with `node_geometry_present 1 ->
  -geometry_point 3`, byte-identical semantic to the Professional
  `.dbs`.
- Blast radius: the fill is a no-op unless `print_node_geometry_present
  -yes` is entered (single corpus test); `geometry()` gained only
  GET_IF_EXISTS lookups of the two new per-geometry records.

## Pending

- Deletion on an empty step fills the index with the LONG_MIN sentinel
  of `db_delete_index` (not a real deallocation); the corpus never
  reads an emptied node record.
