# control_mesh_change_element_group

## Files and functions

- `delete.cc` — `void mesh_delete_keep( long int icontrol )` (line 264). Changes the group of all elements of `group_from` to `group_to`.
- `top.cc` — called in `step_start` (lines 525–530).
- `tochnog.h` — enum `CONTROL_MESH_CHANGE_ELEMENT_GROUP` (line 211).
- `database.cc` — keyword registration (lines 492–495).

## Implementation details

- Reads `long int ceg[2]` with
  `db( CONTROL_MESH_CHANGE_ELEMENT_GROUP, icontrol, ceg, ddum, length, VERSION_NORMAL, GET )`;
  `eg_from = ceg[0]` and `eg_to = ceg[1]`.
- Iterates over all active elements; the group is read with
  `db( ELEMENT_GROUP, ielem, &eg, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS )`
  (`eg` defaults to 0). When `eg == eg_from`, the new group is written with
  `db( ELEMENT_GROUP, ielem, &eg_to, ddum, ldum, VERSION_NORMAL, PUT )`
  (length 1).
- Closes with `mesh_has_changed( VERSION_NORMAL )` although no element is deleted.
- `database.cc`: `type = INTEGER`, `data_length = 2`, `data_class = CONTROL`.

## External dependencies

None. Uses the internal database API (`db`, `db_active_index`, `db_max_index`).

## Hardcoded parameters / pending refactorings

- The change-group record has a fixed length of 2 (from, to) while the delete/keep records use `DATA_ITEM_SIZE`.
- The routine is invoked through the generic delete/keep path although it deletes nothing; a dedicated `mesh_change_element_group` would be cleaner.
- Elements without `ELEMENT_GROUP` are treated as group 0.
