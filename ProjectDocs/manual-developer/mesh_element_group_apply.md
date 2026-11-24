# mesh_element_group_apply — developer

- Records `MESH_ELEMENT_GROUP_APPLY` (INTEGER, up to DATA_ITEM_SIZE,
  record index) and `CONTROL_MESH_ELEMENT_GROUP_APPLY` (INTEGER, CONTROL),
  `database.cc` after `ELEMENT_GROUP`.
- `element_group_applied( group )` (`group.cc`): control record of the
  current `ICONTROL` first, else the union of the global records, else 1.
- `elem()` (`elem.cc`) returns right after the `group_time` test when the
  group is not applied; `bounda_node_dormant()` (`bounda.cc`) marks nodes
  of only non-applied (or not alive) groups dormant; its early return now
  also checks `element_group_apply_active()`.
