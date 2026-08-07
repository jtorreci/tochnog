# control_mesh_remove

## Files and functions

- `delete.cc` — `void mesh_remove( long int icontrol )` (line 354). Reads the keyword record, selects the method and deletes the matching elements.
- `delete.cc` — `void delete_element( long int element, long int version )` (line 26). Deletes every `ELEMENT`-class item of the given element via `db_delete_index`.
- `top.cc` — called in `step_start` (lines 550–553) when `db_active_index( CONTROL_MESH_REMOVE, icontrol, VERSION_NORMAL )` is active.
- `tochnog.h` / `tochnog-mod.h` — enum `CONTROL_MESH_REMOVE` (lines 258 / 251), methods `METHOD1` and `METHOD2` (lines 732–733 / 725–726), prototype `void mesh_remove( long int icontrol )` (line 1371).
- `database.cc` — keyword registration (lines 760–764).

## Implementation details

- Reads the record with `db( CONTROL_MESH_REMOVE, icontrol, list, ddum, length, VERSION_NORMAL, GET )`; `list[0]` is the method, `list[1..]` the element groups, with `ngrp = length - 1`. Returns early if `length < 1` or `max_elem < 0`.
- `-method1` (`method == -METHOD1`): iterates all active elements; for each element whose `ELEMENT_GROUP` equals `list[1]`, checks whether every node of it is also a node of some element belonging to the groups `list[2..]` (node connectivity comparison). If so, removes it with `delete_element`. Requires `ngrp >= 1`.
- Any other method prints `"Error: control_mesh_remove method must be -method1 or -method3."` and exits with `TN_EXIT_STATUS`.
- Closes with `mesh_has_changed( VERSION_NORMAL )`.
- `database.cc`: `type = INTEGER`, `data_length = DATA_ITEM_SIZE`, `fixed_length = 0`, `data_class = CONTROL`.

## External dependencies

None. Uses the internal database API (`db`, `db_active_index`, `db_max_index`) and `delete_element`.

## Hardcoded parameters / pending refactorings

- `-method3` (remove elements where all nodes have an mpc) is not implemented: `NODE_MPC` does not exist and `METHOD3` is not in the enum — only `METHOD1`/`METHOD2` are defined.
- `-method1` uses node-coincidence comparison, not a real geometric containment test.
- `METHOD2` is declared but unused by `mesh_remove`.
- `list[DATA_ITEM_SIZE]` limits the number of groups usable in a single record.
