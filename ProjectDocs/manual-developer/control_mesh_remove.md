# control_mesh_remove

## Files and functions

- `mesh_remove.cc` — `void mesh_remove( long int icontrol )` (moved from `delete.cc` in F5 T4 U9e). Reads the keyword record, selects the method and deletes (or, with `control_mesh_remove_really -no`, deactivates) the matching elements; see [control_mesh_remove_really](control_mesh_remove_really.md).
- `delete.cc` — `void delete_element( long int element, long int version )` (line 26). Deletes every `ELEMENT`-class item of the given element via `db_delete_index`.
- `top.cc` — called in `step_start` (every step) when `db_active_index( CONTROL_MESH_REMOVE, icontrol, VERSION_NORMAL )` is active.
- `tochnog.h` / `tochnog-mod.h` — enum `CONTROL_MESH_REMOVE`, methods `METHOD1`, `METHOD2` and `METHOD3` (F5 T4 U9f, the word `-method3` registered in `database.cc`), prototype `void mesh_remove( long int icontrol )`.
- `point_el.cc` — `point_el()` (the `post_point` locator), reused unchanged for the `-method1` containment test.
- `database.cc` — keyword registration (lines 760–764).

## Implementation details

- Reads the record with `db( CONTROL_MESH_REMOVE, icontrol, list, ddum, length, VERSION_NORMAL, GET )`; `list[0]` is the method, `list[1..]` the element groups, with `ngrp = length - 1`. Returns early if `length < 1` or `max_elem < 0`.
- `-method1` (`method == -METHOD1`, F5 T4 U9f): the reference elements are collected once per pass (`mesh_remove_reference_add`): active elements of the groups `list[2..]` except `list[1]` (Professional: elements of `element_group_0` are no containers), with their `NODE` coordinates, an eps_iso (2.02e-4 for tensor-product elements, 0.85e-4 for triangle / tetrahedron / prism families, values inside the measured Professional windows) and a bounding box widened by 1 % of its size as a prefilter. Every element (active, and the inactive ones kept in `ELEMENT_REMOVE_INACTIVE`, with their kept coordinates when a node no longer exists) of group `list[1]` is removed when each node lies inside some reference element other than itself (`mesh_remove_inside_reference`, `point_el( point, coords, weight, name, nnol, eps_iso )` on the full node list, so quadratic references are curved). Requires `ngrp >= 1`; the element is deleted (`_really -yes`, default) or deactivated (`-no`).
- `-method3` (`method == -METHOD3`): a node flag array from the active `MPC_NODE_NUMBER` records (`mpc[0]`, the constrained node; the same rule as `geometry_node_has_mpc` in geometry.cc, a static helper there); an element is removed when all its nodes are flagged. Groups after `-method3` stop with `db_error( CONTROL_MESH_REMOVE, icontrol )` (= Professional).
- Closes with `mesh_has_changed( VERSION_NORMAL )` (after every `-method1` pass with groups, as before, or when an element was deleted, deactivated or activated; a `-method3` pass only when something changed).
- `database.cc`: `type = INTEGER`, `data_length = DATA_ITEM_SIZE`, `fixed_length = 0`, `data_class = CONTROL`.

## control_mesh_remove_geometry (F5 T4 U5c)

`mesh_remove()` reads `CONTROL_MESH_REMOVE_GEOMETRY` (INTEGER, 2: the
geometry entity and index, requires `CONTROL_MESH_REMOVE`) and skips a
candidate element unless `element_in_geometry( ielem, entity, -ALL,
NODE_START_REFINED, PROJECT_EXACT, NULL, VERSION_NORMAL )` (geometry.cc,
the shared helper of F5 T4 U4) holds. The reference elements of the
other groups are not filtered. Verified on the Professional 02-08-2026
(the two geometries of `remove_geometry1` / `remove_geometry2`).

## External dependencies

None. Uses the internal database API (`db`, `db_active_index`, `db_max_index`), `delete_element` and `point_el`.

## Differences with the Professional version

- Curved quadratic reference elements: the GNU tests the curved element, the Professional its straight corner element (DIVERGENCES REMOVE-METHOD1-CURVED, `control_mesh_remove_method1_curved1`).

## Hardcoded parameters / pending refactorings

- The `-method1` test is O(candidate nodes x reference elements), with the bounding-box prefilter; no spatial search tree.
- The containment uses `NODE` (the current coordinates), the `_geometry` restriction `NODE_START_REFINED`.
- `-method3` reads the `mpc_node_number` records present when `step_start` runs; the generated ones are made by `mpc_node_apply` during the iterations of a step (top.cc, bounda.cc), so from the code a remove on the first step sees only the input records (not tested).
- `METHOD2` is declared but unused by `mesh_remove`.
- `list[DATA_ITEM_SIZE]` limits the number of groups usable in a single record.
