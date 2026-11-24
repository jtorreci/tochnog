# control_zip — developer

- Record `CONTROL_ZIP` (INTEGER, CONTROL; `database.cc` after `ZIP`).
- `top()` (`top.cc`): at the end of each control index (after
  `control_system_call`), `-YES` calls `zip_output_files()` (`miscel.cc`,
  the shell loop factored out of the global `zip` of `exit_tn()`); a
  non-zero status gives a warning.
