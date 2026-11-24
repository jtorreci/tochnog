# group_materi_elasti_young_user (developer)

## Where

- `database.cc`: `GROUP_MATERI_ELASTI_YOUNG_USER` (INTEGER, 1, MATERI);
  enums in `tochnog.h` / `tochnog-mod.h`.
- `stress.cc` `set_stress()`: the `group_materi_elasti_young` block also
  runs when the switch is -yes; `user_young( new_unknowns, old_hisv,
  rotated_old_sig, young )` sets the modulus before `C_matrix`.
- `user.cc` `user_young()`: stub, error and exit; prototype in
  `tochnog.h`.

## Pending

- The young_polynomial / young_power laws do not call the routine.
