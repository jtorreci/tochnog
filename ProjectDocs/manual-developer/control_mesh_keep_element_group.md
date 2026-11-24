# control_mesh_keep_element_group

## Files and functions

- `delete.cc` — `void mesh_delete_keep( long int icontrol )` (line 264). Keeps the elements of the listed groups and deletes the rest.
- `delete.cc` — `void delete_element( long int element, long int version )` (line 26).
- `math.cc` — `long int array_member( long int list[], long int i, long int n, long int &indx )` (line 52).
- `top.cc` — called in `step_start` (lines 525–530).
- `tochnog.h` — enum `CONTROL_MESH_KEEP_ELEMENT_GROUP` (line 232).
- `database.cc` — keyword registration (lines 611–615).

## Implementation details

- Reads the group list via `db( CONTROL_MESH_KEEP_ELEMENT_GROUP, icontrol, list, ddum, nlist, VERSION_NORMAL, GET )`.
- Iterates over all active elements; the element group is read with
  `db( ELEMENT_GROUP, ielem, &eg, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS )`,
  where `eg` is initialized to 0, so elements without an explicit `ELEMENT_GROUP`
  are treated as group 0.
- Elements whose group is not in the list are removed with `delete_element`.
- Closes with `mesh_has_changed( VERSION_NORMAL )`.
- `database.cc`: `type = INTEGER`, `data_length = DATA_ITEM_SIZE`, `fixed_length = 0`, `data_class = CONTROL`.

## External dependencies

None. Uses the internal database API (`db`, `db_active_index`, `db_max_index`), `array_member` and `delete_element`.

## Hardcoded parameters / pending refactorings

- Elements without `ELEMENT_GROUP` default to group 0 and are therefore removed unless 0 is in the keep list.
- The group is fetched with `GET_IF_EXISTS`; the code relies on `eg` being initialized to 0 before the read.
