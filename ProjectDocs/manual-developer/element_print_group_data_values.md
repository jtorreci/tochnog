# print_group_data + element_print_group_data_values (developer)

## Where

- `group.cc`, `element_print_group_data_values_calculate()`: called from
  `top.cc` `step_start()` (serial, after the geometry-present scans). For
  each element: its `ELEMENT_GROUP` (default 0), then for every name of
  `PRINT_GROUP_DATA` whose record name starts with `group_` the values of
  that record at the group index (INTEGER records converted to double); an
  absent record gives one 0. Real records go through `get_group_data()`
  (F5 T4 U3, see dependency_item.md): unknowns = element node average,
  `-time_current` = `time_current` argument (end of the starting step).
- `database.cc`: `PRINT_GROUP_DATA` widened from length 1 to
  `DATA_ITEM_SIZE` (variable); `ELEMENT_PRINT_GROUP_DATA_VALUES`
  DOUBLE_PRECISION, `DATA_ITEM_SIZE`, variable, `print_only`, class ELEMENT;
  the name `element_print_group_data_value` (Professional spelling) maps to
  it in the alias table of `db_number()`.

## Pending

- No GiD output of the values (manual 6.991 mentions it).
- Test `print_group_data1`.
