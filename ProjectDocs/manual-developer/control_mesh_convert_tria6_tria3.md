# control_mesh_convert_tria6_tria3

## Files and functions

- `mesh.cc` — `mesh_convert_tria6_tria3()`: called from `top.cc`
  `step_start()` on every step, after `mesh_convert_hex20()`.
- `database.cc` — `CONTROL_MESH_CONVERT_TRIA6_TRIA3` (INTEGER switch,
  CONTROL; F5 T4 U5b). Enum after `CONTROL_MESH_CONVERT_ELEMENT_GROUP` in
  `tochnog.h` / `tochnog-mod.h`.

## Implementation details

When `ICONTROL` has the record with `-yes`, every `-tria6` (length 7) is
rewritten in place as `-tria3 el[1] el[3] el[6]` (corners n0, n2, n5 of
the tria6 order); then the nodes no element references are deleted
(`delete_node`) and `mesh_has_changed()` is called. Idempotent: later
steps find no tria6. Element-class records of the element are kept.

Verified against the Professional 02-08-2026 (probe on two tria6: same
elements and nodes).

## Pending

- The element-class data of a converted element (for example an
  `element_dof` record written before the conversion) is not resized;
  the record is meant for the mesh set-up, before the calculation.
