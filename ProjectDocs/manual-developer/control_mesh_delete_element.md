# control_mesh_delete_element

## Files and functions

- `delete.cc` — `void mesh_delete_keep( long int icontrol )` (line 264). Reads the element list and calls `delete_element` for each active listed element.
- `delete.cc` — `void delete_element( long int element, long int version )` (line 26). Deletes every `ELEMENT`-class item of the given element via `db_delete_index`.
- `top.cc` — called in `step_start` (lines 525–530) when the keyword is active.
- `tochnog.h` — enum `CONTROL_MESH_DELETE_ELEMENT` (line 217) and prototype (line 1365).
- `database.cc` — keyword registration (lines 529–533).

## Implementation details

- Part of `mesh_delete_keep`, shared with keep_element, keep_element_group and change_element_group; the branch is selected by `db_active_index`.
- Reads the integer list with `db( CONTROL_MESH_DELETE_ELEMENT, icontrol, list, ddum, nlist, VERSION_NORMAL, GET )`; `nlist` is filled with the number of items in the record.
- For each entry, `db_active_index( ELEMENT, list[i], VERSION_NORMAL )` guards against non-existing element numbers.
- Closes with `mesh_has_changed( VERSION_NORMAL )`.
- `database.cc`: `type = INTEGER`, `data_length = DATA_ITEM_SIZE`, `fixed_length = 0`, `data_class = CONTROL`.

## External dependencies

None. Uses the internal database API (`db`, `db_active_index`) and `delete_element`.

## Hardcoded parameters / pending refactorings

- The four delete/keep/group keywords share one routine; the behavior is switched only by `db_active_index` checks. Could be split into dedicated functions.
- `list[DATA_ITEM_SIZE]` limits the number of elements deletable in a single record.
