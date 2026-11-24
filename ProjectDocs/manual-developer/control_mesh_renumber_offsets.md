# control_mesh_renumber_element_group_offset / _element_geometry_offset

## Files and functions

- `top.cc` — the `CONTROL_MESH_RENUMBER` block (after `renumbering()`):
  reads `CONTROL_MESH_RENUMBER_ELEMENT_GROUP_OFFSET` and
  `CONTROL_MESH_RENUMBER_ELEMENT_GEOMETRY_OFFSET` and adds them to
  `ELEMENT_GROUP` / `ELEMENT_GEOMETRY` of every active element (a missing
  record is read as 0 and then written).
- `database.cc` — both INTEGER, 1, CONTROL, `data_required =
  CONTROL_MESH_RENUMBER` (F5 T4 U5c).

## Implementation details

The offsets are applied only by `control_mesh_renumber`; the internal
`renumbering()` calls (`control_mesh_split`, refinement) do not see them.
The records are read with explicit names (the buffer audit cannot follow
a record number held in a variable). Verified on the Professional
02-08-2026 (groups 0 1 2 -> 100 101 102, element_geometry 0 3 0 -> 7 10 7).
