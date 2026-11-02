# control_check_data — developer

- Record `CONTROL_CHECK_DATA` (INTEGER, CONTROL; `database.cc` after
  `CHECK_DATA`).
- `top()` (`top.cc`): at the entry of each control index (after
  `ICONTROL` is set) `-YES` calls `check_data_integrity()`
  (`database.cc`), the routine of the global `check_data`.
