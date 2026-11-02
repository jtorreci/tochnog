# control_mesh_copy

## Files and functions

- `mesh.cc` — `void mesh_copy( double move_coords[] )` (line 289). Duplicates the mesh displaced over `move_coords`.
- `top.cc` — called in `step_start` (lines 516–522) when `db_active_index( CONTROL_MESH_COPY, icontrol, VERSION_NORMAL )` is true.
- `tochnog.h` — enum `CONTROL_MESH_COPY` (line 217) and prototype (line 1358).
- `database.cc` — keyword registration (lines 529–533).

## Implementation details

- Same duplication pattern as `mesh_mirror`, but instead of negating the mirror
  coordinate each copied node is moved over the displacement:
  `new_coords[i] = coords[i] + move_coords[i]` for `i < ndim`.
- Every active node is duplicated with index `n + max_node_old + 1`; the same
  displacement is applied to `NODE_START_REFINED`, and `NODE_DOF` is copied via
  `db_dbl( NODE_DOF, inod, ... )` with length `nuknwn`.
- Elements are duplicated with `create_element( ielem, new_elem, new_nodes, length, ... )`,
  which copies all element items; new element nodes are the old ones offset by
  `max_node_old + 1`.
- Closes with `mesh_has_changed( VERSION_NORMAL )`.
- `database.cc`: `type = DOUBLE_PRECISION`, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0`, `data_class = CONTROL`.

## External dependencies

None. Uses the internal database API (`db`, `db_dbl`, `db_active_index`,
`db_max_index`), `create_element` and `mesh_has_changed`.

## Hardcoded parameters / pending refactorings

- In `top.cc` the copy values are read into `control_copy[MDIM]` and only the
  first `ndim` are used; `copy_length` is read but not passed to `mesh_copy`.
- The node duplication loop is duplicated with `mesh_mirror`; a shared helper
  (e.g. `duplicate_mesh(coord_transform)`) could unify both.

## Once per control index (F4 T7)

The call in `top.cc` `step_start()` is guarded by `first_of_control`
(see [control_mesh_once_per_control](control_mesh_once_per_control.md)).

## Node and element data of the duplicate (F4 T8, BACKLOG C45)

- `mesh.cc` `mesh_duplicate_node_data( inod, new_node )` (static): copies
  every node-class item except `NODE` / `NODE_START_REFINED` (written by
  the caller with the coordinate transformation). Before, only `NODE_DOF`
  was copied; `renumber()` reads `NODE_DOF_START_REFINED` of every node and
  stopped with "Error detected for data item : node_dof_start_refined".
  The source is copied into a local buffer before the PUT: a PUT beyond
  the allocated range reallocates the data array, so a `db_dbl()` pointer
  passed to it dangled (latent use-after-free, also in `mesh_rotate_3d`,
  which now uses the same helper).
- `mesh.cc` `mesh_duplicate_element_state( ielem, new_elem )` (static):
  copies `ELEMENT_INTPNT_DOF` and `ELEMENT_DOF_INITIALISED` (versions NORMAL and
  NEW) after `create_element()`, which resets them (right for refinement)
  and so re-initialised the duplicate from the smoothed nodal values (1 %
  off the Professional on a loaded block).
- New element numbers: `ielem + max_element + 1` (Professional probes),
  before sequential from `max_element + 1`.
- Verification: probes on a loaded 3x2-node block (gravity, step on index
  0, copy on index 1 with two steps): all nodal dofs = Professional to 4e-10
  (round-off level of the untouched block).

## Buffered node-item copies in mesh.cc (F4 T8, ASan)

`mesh.cc` `mesh_copy_node_item( idat, from, to, len )` (static) copies one
node-class item through a local buffer. It replaces every
`db_dbl()/db_int()` pointer passed straight to a `db(... PUT)` of the
same item in `mesh.cc` (`mesh_duplicate_node_data` - copy, mirror,
rotate_3d - and `mesh_convert_quad8`, `mesh_convert_hex20`,
`mesh_extrude` x3): a PUT beyond the allocated range reallocates the
array and the source pointer then reads freed memory (heap-use-after-
free reported by the F4 T9 ASan run). Results unchanged (all 297 suite
`.dbs` identical apart from timing lines; ASan clean on the copy/mirror
tests and the Professional `quad8`, `hex20`). The same pattern remains in
`interface.cc` `interface_convert()` (lines ~1450 / ~1631; ASan on the
Professional `interface_bar2_hex8`, `interface_bar3_quad8`,
`interface_quad8_hex20`), outside F4 T8. Fixed in F4 T9 (2026-10-10, `interface_copy_node_dbl`, see
`group_interface.md`).

