# print_group_data

## Implementation

- Keyword `PRINT_GROUP_DATA` (INTEGER, no_index, data_length 1,
  fixed_length 1) registered in `database.cc`; enum appended in
  `tochnog.h`/`tochnog-mod.h` (same enum, same order); later widened to
  `DATA_ITEM_SIZE`. Consumed by `element_print_group_data_values_calculate()`
  (see element_print_group_data_values.md).

## Pending

- The GiD group-data plot.
