# control_mesh_keep_node

## Files and functions

- `delete.cc` — `void mesh_delete_keep( long int icontrol )` (line 265). Keeps the listed nodes and deletes the rest.
- `delete.cc` — `void delete_node( long int inod, long int version )` (line 44). Deletes a node and all its node-class items.
- `math.cc` — `long int array_member( long int list[], long int i, long int n, long int &indx )` (line 52). Returns true if `i` is in `list[0..n-1]`.
- `top.cc` — called in `step_start` (lines 534–540) when `db_active_index( CONTROL_MESH_KEEP_NODE, icontrol, VERSION_NORMAL )` is true.
- `tochnog.h` — enum `CONTROL_MESH_KEEP_NODE` (line 234) and `delete_node` prototype (line 1177).
- `database.cc` — keyword registration (lines 623–626).

## Implementation details

- Reads the keep list via `db( CONTROL_MESH_KEEP_NODE, icontrol, list, ddum, nlist, VERSION_NORMAL, GET )`.
- Iterates over all active nodes `0..max_node` (`db_max_index`) and uses `array_member( list, inod, nlist, ldum )` to decide whether the node is kept; non-listed nodes are removed with `delete_node`.
- `delete_node` loops over every data item (`0..MDAT`); for items whose `db_data_class` is `NODE` it calls `db_delete_index( idat, inod, version )`, so the node and all its associated items are removed.
- Closes with `mesh_has_changed( VERSION_NORMAL )`.
- `database.cc`: `type = INTEGER`, `data_length = DATA_ITEM_SIZE`, `fixed_length = 0`, `data_class = CONTROL`.

## External dependencies

None. Uses the internal database API (`db`, `db_active_index`, `db_max_index`, `db_data_class`, `db_delete_index`), `array_member` and `delete_node`.

## Hardcoded parameters / pending refactorings

- No check is done on elements that reference deleted nodes; they keep dangling references and will fail if used (intended for visualization).
- `array_member` reports the found index through a reference (`&indx`) that is discarded here (`ldum`).
- `delete_node` is declared twice in `tochnog.h` (lines 1177 and 1179); one duplicate declaration could be removed.
- Same shared-routine remark as `control_mesh_keep_element`.
