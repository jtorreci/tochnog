# group_time

## Implementation

- `group.cc` `group_time_alive(element_group, time_total)`: the single
  activity test of a `group_time` record: active when there is no record,
  when `time_total == 0 && birth == 0`, or when
  `birth < time_total <= death`. Declared in `tochnog.h` /
  `tochnog-mod.h`.
- `elem.cc` `elem()`: returns early (no element contribution) when the
  group has a `group_time` record and `group_time_alive()` is false. Same
  test as before 2026-10-06, now through the helper.
- `bounda.cc` `bounda_node_dormant(max_node, time_total)` (static,
  C39 2026-10-06): when at least one `group_time` record exists, loops
  over the active elements, marks the nodes of elements whose group is
  alive and returns an array `[0..max_node]` with 1 for a node that has
  elements but none alive ("dormant"). Returns NULL without `group_time`
  records, so models without the record run exactly the old code.
  `bounda()` skips the dormant nodes for the prescribed-dof branch
  (`bounda_dof` / `bounda_unknown`, the `unknown` flag): `NODE_BOUNDED`
  stays 0 and `NODE_DOF` new is not overwritten. `bounda_force` records
  are not affected. The array is freed next to `nodes_in_geometry`.

## Verification

- Professional 02-08-2026, test `other/truss8` (truss born at t = 1,
  node 1 velx 0.1 from t = 0): Pro node 1 disx 0 at t = 0.5/1, 0.05 at
  1.5, 0.1 at 2; GNU before 0.05/0.1/0.15/0.2, after = Pro at every step
  (node 2 = Pro before and after).
- Own test `validation-suite/test-2026/group_time_bounda_dormant1`
  (registered `run_ok` in `scripts/build_safe.sh`): adds a second truss
  without `group_time` whose bounded node keeps moving (disx 0.2). RED on
  the old binary (rc 1), GREEN after (rc 0); the Professional passes it.
- Corpus: only `other/truss8` and the new test contain `group_time`.

## Hardcoded data / pending

- Only `group_time` makes an element inactive here; elements switched
  off by `control_mesh_activate_gravity` method 1 before their birth, or
  removed by `control_mesh_delete*`, are not considered dormant (not
  measured on the Professional).
- Death (BACKLOG C39): the Professional zeroes the velocity of the nodes
  of a dead group (also of the free node 2 of truss8-like probes with
  `group_time 0 0. 1.`); the GNU keeps the old velocity of nodes without
  an equation. Not implemented.
- Nodes without any element are not marked dormant (the Professional
  aborts such a model with a singular-matrix message).
