# control_mesh_duplicate_element_group

## Files and functions

- `delete.cc` — `mesh_duplicate_element_group(icontrol)`; called from
  `top.cc` `step_start()` after `mesh_delete_keep()`, when the running
  control index has the record.
- `database.cc` — `CONTROL_MESH_DUPLICATE_ELEMENT_GROUP` (INTEGER, 2,
  CONTROL; F5 T4 U5c). Enum after `CONTROL_MESH_REMOVE_GEOMETRY` in
  `tochnog.h` / `tochnog-mod.h`.

## Implementation details

For every active element up to the highest element number whose
`ELEMENT_GROUP` (default 0) equals the old group, `create_element( ielem,
++new_elem, el, ... )` copies the element record and every element-class
record, then `ELEMENT_GROUP` of the copy is set to the new group;
`mesh_has_changed()` when anything was copied. Verified on the
Professional 02-08-2026 (copy numbering, groups and the copied
`element_geometry`).

## Pending

- Each pass of the control index makes new copies (a record at an index
  with `control_timestep` duplicates again at every step, as the other
  per-index mesh commands of `step_start`).
