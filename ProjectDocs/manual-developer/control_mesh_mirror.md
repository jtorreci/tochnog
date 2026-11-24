# control_mesh_mirror

## Files and functions

- `mesh.cc` — `void mesh_mirror( long int axis )` (line 229). Duplicates the mesh mirrored about the x=0, y=0 or z=0 plane.
- `top.cc` — called in `step_start` (lines 509–514) when `db_active_index( CONTROL_MESH_MIRROR, icontrol, VERSION_NORMAL )` is true.
- `tochnog.h` — enum `CONTROL_MESH_MIRROR` (line 237) and prototype (line 1357).
- `database.cc` — keyword registration (lines 639–642).

## Implementation details

- `axis` is validated against the program constants `-X`/`-Y`/`-Z` and mapped to
  `axis_i` 0/1/2; any other value is a fatal error (`exit(TN_EXIT_STATUS)`).
- Every active node is duplicated with `new_coords[i] = coords[i]` and
  `new_coords[axis_i] = -coords[axis_i]`. The copy of node `n` gets index
  `n + max_node_old + 1`, where `max_node_old` is the max node index before duplication.
- `NODE_START_REFINED` is copied with the same mirrored coordinate when active.
- `NODE_DOF` is copied via `db_dbl( NODE_DOF, inod, ... )` with length `nuknwn`.
- Elements are duplicated with `create_element( ielem, new_elem, new_nodes, length, ... )`,
  which copies ALL element items (`element_dof`, `element_dof_initialised`,
  `nonlocal`, etc.). New element nodes are the old ones offset by `max_node_old + 1`.
- Closes with `mesh_has_changed( VERSION_NORMAL )`.
- `database.cc`: `type = INTEGER`, `data_length = 1`, `data_class = CONTROL`.

## External dependencies

None. Uses the internal database API (`db`, `db_dbl`, `db_active_index`,
`db_max_index`), `create_element` and `mesh_has_changed`.

## Hardcoded parameters / pending refactorings

- Only the first `ndim` coordinate components are mirrored, and only about the
  plane at 0 (no offset mirror plane).
- The node/`NODE_START_REFINED`/`NODE_DOF` duplication loop is duplicated with
  `mesh_copy`; a shared duplication helper could be extracted.

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
  0, mirror on index 1 with two steps): all nodal dofs = Professional to 4e-10
  (round-off level of the untouched block).
