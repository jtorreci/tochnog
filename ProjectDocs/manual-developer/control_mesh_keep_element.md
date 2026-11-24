# control_mesh_keep_element

## Files and functions

- `delete.cc` — `void mesh_delete_keep( long int icontrol )` (line 264). Keeps the listed elements and deletes the rest.
- `delete.cc` — `void delete_element( long int element, long int version )` (line 26).
- `math.cc` — `long int array_member( long int list[], long int i, long int n, long int &indx )` (line 52). Returns true if `i` is in `list[0..n-1]`.
- `top.cc` — called in `step_start` (lines 525–530).
- `tochnog.h` — enum `CONTROL_MESH_KEEP_ELEMENT` (line 231) and prototype (line 1365).
- `database.cc` — keyword registration (lines 605–609).

## Implementation details

- Reads the keep list via `db( CONTROL_MESH_KEEP_ELEMENT, icontrol, list, ddum, nlist, VERSION_NORMAL, GET )`.
- Iterates over all active elements `0..max_elem` (`db_max_index`) and uses `array_member( list, ielem, nlist, ldum )` to decide whether the element is kept; non-listed elements are removed with `delete_element`.
- Closes with `mesh_has_changed( VERSION_NORMAL )`.
- `database.cc`: `type = INTEGER`, `data_length = DATA_ITEM_SIZE`, `fixed_length = 0`, `data_class = CONTROL`.

## External dependencies

None. Uses the internal database API (`db`, `db_active_index`, `db_max_index`), `array_member` and `delete_element`.

## Hardcoded parameters / pending refactorings

- `array_member` reports the found index through a reference (`&indx`) that is discarded here (`ldum`).
- Same shared-routine remark as `control_mesh_delete_element`.
