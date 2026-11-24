# control_data_arithmetic (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enums `CONTROL_DATA_ARITHMETIC`,
  `CONTROL_DATA_ARITHMETIC_DOUBLE`; switch enums `PLUS`, `MINUS`,
  `DIVIDE` added next to the existing `MULTIPLY` (headers in sync).
- `database.cc` — registrations (arithmetic: `INTEGER`,
  `DATA_ITEM_SIZE`, fixed_length 0; double: `DOUBLE_PRECISION`, 1,
  cross-required) and the switch name registrations `plus`/`minus`/
  `divide` next to `multiply` (switches resolve through the same
  `db_number` name table as every keyword).
- `data.cc` — new block in `data()`, fires at each `step_close` whose
  `ICONTROL` matches the record index. Record layout
  `[ item, index | -RA range..., number | -ALL, operat ]`: number and
  operat are the LAST two slots; the index is a single value or a
  `-ra ... -ra` range expanded with `range_expand` (same pattern as
  `control_data_delete`). `-ALL` loops over all values of the record.
  Integer-typed items are rejected (consistent with `change_dataitem`);
  division by zero is rejected up front. A negative `number` is resolved
  through `DOF_LABEL` like `change_dataitem` does (own local buffer —
  the shared `dof_label` pointer is not allocated yet at that point).

## Timing (gotcha)

`data()` runs in `step_close` at the END of each timestep: the change is
applied once per STEP of the owning `control_timestep` block. And note
`control_timestep i dt increment`: the SECOND number is the increment
DURATION, not the end time — `0.1 0.2` from t=0.1 is TWO steps (to
0.3), so a multiply record fires twice. Verified empirically (1000 ->
4000 instead of 2000); the test uses one-step blocks.

## Verification

Test `cda_arith` (suite 67/67): young 1000; block 1 `-multiply 2.0` on
number 0; block 2 `-plus 500.0` with `-all`; final record value read by
a generic `target_item` = 2500 exactly. Covers single index, plain
number, `-all`, multiply and plus paths.
