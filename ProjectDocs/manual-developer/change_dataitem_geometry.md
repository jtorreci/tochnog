# change_dataitem_geometry (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — data item enum
  `CHANGE_DATAITEM_GEOMETRY`.
- `database.cc` — registration: `INTEGER`, `data_length = 2`
  (geometry_entity_name, geometry_entity_index),
  `data_required = CHANGE_DATAITEM`.
- `data.cc` — in the `change_dataitem` loop, when the record exists and
  the target is a `group_*` item: group-split materialization.

## Design: group split

A `group_*` record is one database record per group, so "apply the
change only to some elements of the group" cannot be done by editing the
record. Instead the group is SPLIT:

1. On the first application (found=1), a clone group index is allocated
   (max group index + 1) and remembered in a function-static map
   `clone_map[original_group] -> clone_group` (grow-on-demand array).
2. Every data item whose name starts with `group_` and that is active at
   the original group index is copied to the clone (integer and double
   variants; `group_type` itself included — it also starts with
   `group_`). Buffer safety: all group records have
   `data_length <= DATA_ITEM_SIZE`.
3. Nodes inside the geometry are found with `geometry()` (same call
   shape as the `CONTROL_DATA_INITELDOF_GEOMETRY` block; coordinates
   from `NODE_START_REFINED`). Elements of the original group with ALL
   nodes inside get their `ELEMENT_GROUP` record rewritten to the clone
   (same PUT pattern as delete.cc's group reassignment).
4. `data_item_index` is redirected to the clone, so the value change
   (this and every later timestep) applies to the clone record only.

Element state (NODE_DOF stresses/history) is untouched: moving an element
between groups changes which parameters it reads, not its state — the
same semantics as `control_mesh_change_element_group`.

## Limitations (documented)

- Elements must be FULLY inside the geometry; partially covered elements
  keep the original parameters (conservative zoning).
- `clone_map` is process state: after a restart the mapping is lost and a
  continued run would split again (creating a second clone). Restarting
  mid-change with geometry restriction is not supported.
- Only ONE geometry-restricted change per original group is mapped; a
  second change record targeting the same group with a different
  geometry shares the first clone (first restriction wins for element
  placement; both apply values to the same clone).
- Region selectors: `geometry_brick` (inside-a-box) is the natural
  choice. Note `geometry_circle` means ON the circle for non-delete
  projection types (only CONTROL_MESH_DELETE_GEOMETRY gets the
  inside-the-circle semantics) — see geometry.cc.

## Verification

Test `cd_geom` (suite 64/64): two disjoint quad4 blocks in pure shear,
both in group 0, capped by mohr_coul_direct (c=1 -> sigxy ~ 1 both).
`geometry_brick 1` covers only the left block; from t=0.05 c -> 0 with
the restriction. After the split: left block sigxy ~ 0 (clone c=0),
right block sigxy ~ 1 (original group), and a generic target reads the
original record still holding c=1.0 — the restriction is proven, not
just the absence of a crash.
