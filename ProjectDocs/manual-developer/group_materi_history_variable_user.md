# group_materi_history_variable_user / _parameters (developer)

## Where

- `database.cc`: `GROUP_MATERI_HISTORY_VARIABLE_USER` (INTEGER, 1) and
  `GROUP_MATERI_HISTORY_VARIABLE_USER_PARAMETERS` (DOUBLE list), MATERI.
- `materi.cc` `materi()`: after `set_stress()` and `materi_stress_null()`,
  `user_materi_history_variable( parameters, nparameters, new_unknowns,
  old_hisv, new_hisv, new_sig )` when the switch is -yes.
- `user.cc`: stub, error and exit; prototype in `tochnog.h`.
