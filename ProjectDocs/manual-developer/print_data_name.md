# print_data_name — developer

- Record `PRINT_DATA_NAME` (INTEGER, `no_index`; `database.cc` after
  `MESSAGE`).
- `top()` (`top.cc`), after the introduction message: writes
  `db_name(idat)` for `idat` 0..MDAT-1 to `tochnog_data_name.txt`,
  skipping empty and blank names (Professional layout: "\nNames:\n" and
  one name per line, enum order, not sorted).
