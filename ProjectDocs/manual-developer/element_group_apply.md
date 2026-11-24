# element_group_apply — developer

- Records `ELEMENT_GROUP_APPLY` (INTEGER, up to DATA_ITEM_SIZE, data
  class ELEMENT), `CONTROL_ELEMENT_GROUP_APPLY` and
  `CONTROL_ELEMENT_GROUP` (INTEGER, CONTROL).
- `control_element_group_apply( icontrol )` (`group.cc`), called by
  `top()` at the entry of each control index (after
  `control_check_data`): PUTs `ELEMENT_GROUP` from position `number`
  (`db_error` when missing), then `area_element_group( VERSION_NORMAL )`
  for `control_element_group -yes`.
