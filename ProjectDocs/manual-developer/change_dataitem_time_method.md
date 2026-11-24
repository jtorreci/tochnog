# change_dataitem_time_method (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — switch enums `COSINUS`, `SINUS`,
  `TANGENT` (alphabetical positions; headers in sync) and data item enum
  `CHANGE_DATAITEM_TIME_METHOD`.
- `database.cc` — keyword registration (`INTEGER`, `data_length = 1`,
  `data_required = CHANGE_DATAITEM`) and the switch NAME registrations
  (`cosinus`, `sinus`, `tangent`) — the parser resolves `-tangent`
  through the same name table as every other switch (db_number).
- `data.cc` — in the `change_dataitem` loop, after the time value `val`
  is computed and before it is applied:
  `val = acos/asin/atan(val)` when the method record exists.

## Semantics

The time table holds trigonometric VALUES; the stored parameter is the
inverse angle. `atan` is the phi-c reduction case: tables give tan(phi),
the stored parameter is phi in radians.

## Verification

Test `cd_method` (suite 64/64), an exact unit test of the transform:
base model of `materi_direct_mc` (pure shear, cap c=1 -> sigxy ~ 1) with

```
change_dataitem              10 -group_materi_plasti_mohr_coul_direct 0 0 -use
change_dataitem_time         10 0.0 0.57735027 100. 0.57735027
change_dataitem_time_method  10 -tangent
```

and a GENERIC target reading the stored record directly:
`target_item 2 -group_materi_plasti_mohr_coul_direct 0 0` = 0.5235988
(= atan(0.57735027), tolerance 1e-3). sigxy stays ~1 (sig_n ~ 0, the cap
is still c=1), proving no side effects on the physics.
