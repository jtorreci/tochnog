# print_group_data

## Implementation

- Keyword `PRINT_GROUP_DATA` (INTEGER, no_index, data_length 1,
  fixed_length 1) registered in `database.cc`; enum appended in
  `tochnog.h`/`tochnog-mod.h` (same enum, same order). Parse-only.

## Pending

- The GiD group-data writing and the `element_print_group_data` fill
  (would need the element_group -> group_* data lookup per element and
  the GiD layer hooks).
