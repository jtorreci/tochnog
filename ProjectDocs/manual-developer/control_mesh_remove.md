# control_mesh_remove

## Files and functions

- `mesh_remove.cc` — `void mesh_remove( long int icontrol )` (moved from `delete.cc` in F5 T4 U9e). Reads the keyword record, selects the method and deletes (or, with `control_mesh_remove_really -no`, deactivates) the matching elements; see [control_mesh_remove_really](control_mesh_remove_really.md).
- `delete.cc` — `void delete_element( long int element, long int version )` (line 26). Deletes every `ELEMENT`-class item of the given element via `db_delete_index`.
- `top.cc` — called in `step_start` (every step) when `db_active_index( CONTROL_MESH_REMOVE, icontrol, VERSION_NORMAL )` is active.
- `tochnog.h` / `tochnog-mod.h` — enum `CONTROL_MESH_REMOVE` (lines 258 / 251), methods `METHOD1` and `METHOD2` (lines 732–733 / 725–726), prototype `void mesh_remove( long int icontrol )` (line 1371).
- `database.cc` — keyword registration (lines 760–764).

## Implementation details

- Reads the record with `db( CONTROL_MESH_REMOVE, icontrol, list, ddum, length, VERSION_NORMAL, GET )`; `list[0]` is the method, `list[1..]` the element groups, with `ngrp = length - 1`. Returns early if `length < 1` or `max_elem < 0`.
- `-method1` (`method == -METHOD1`): iterates all elements (active, and since U9e the inactive ones kept in `ELEMENT_REMOVE_INACTIVE`); for each element whose group equals `list[1]`, checks whether every node of it is also a node of some active element belonging to the groups `list[2..]` that is not removed earlier in the same pass (node connectivity comparison; the result of the former in-place deletion). If so, removes it: `delete_element` with `control_mesh_remove_really -yes` (default), deactivation with `-no`. Requires `ngrp >= 1`.
- Any other method prints `"Error: control_mesh_remove method must be -method1 or -method3."` and exits with `TN_EXIT_STATUS`.
- Closes with `mesh_has_changed( VERSION_NORMAL )` (after every `-method1` pass with groups, as before, or when an element was deleted, deactivated or activated).
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

None. Uses the internal database API (`db`, `db_active_index`, `db_max_index`) and `delete_element`.

## Hardcoded parameters / pending refactorings

- `-method3` (remove elements where all nodes have an mpc) is not implemented: `NODE_MPC` does not exist and `METHOD3` is not in the enum — only `METHOD1`/`METHOD2` are defined.
- `-method1` uses node-coincidence comparison, not a real geometric containment test.
- `METHOD2` is declared but unused by `mesh_remove`.
- `list[DATA_ITEM_SIZE]` limits the number of groups usable in a single record.
