# print_group_data + element_print_group_data_values (developer)

## Where

- `group.cc`, `element_print_group_data_values_calculate()`: called from
  `top.cc` `step_start()` (serial, after the geometry-present scans). For
  each element: its `ELEMENT_GROUP` (default 0), then for every name of
  `PRINT_GROUP_DATA` whose record name starts with `group_` the values of
  that record at the group index (INTEGER records converted to double); an
  absent record gives one 0.
- `database.cc`: `PRINT_GROUP_DATA` widened from length 1 to
  `DATA_ITEM_SIZE` (variable); `ELEMENT_PRINT_GROUP_DATA_VALUES`
  DOUBLE_PRECISION, `DATA_ITEM_SIZE`, variable, `print_only`, class ELEMENT;
  the name `element_print_group_data_value` (Professional spelling) maps to
  it in the alias table of `db_number()`.

## Pending

- The values are the raw group records, not the `get_group_data()` values
  (element distributions and dependency diagrams are not applied).
- No GiD output of the values (manual 6.991 mentions it).
- Test `print_group_data1`.
