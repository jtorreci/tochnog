# control_data_copy / control_data_copy_index (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enums `CONTROL_DATA_COPY`,
  `CONTROL_DATA_COPY_FACTOR`, `CONTROL_DATA_COPY_INDEX`,
  `CONTROL_DATA_COPY_INDEX_FACTOR`.
- `database.cc` — registrations (copy: `INTEGER`, 2; copy_index:
  `INTEGER`, 4; factors: `DOUBLE_PRECISION`, 1, cross-required).
- `data.cc` — static helper `data_copy_apply( from, index_from, to,
  index_to, factor, control_item, icontrol )` doing the single-record
  copy: integer->integer requires factor 1; double->double multiplies
  every value by the factor; mixed types are an input error. Buffers of
  `DATA_ITEM_SIZE` (all records fit). `control_data_copy` loops over all
  active indices of the source item (`db_max_index`); `copy_index` uses
  the explicit from/to indices.

The manual's special case node_inertia -> node_force with factor -1
(d'alembert) is covered by the generic double copy.

## Verification

Test `cda_copy` (suite 67/67): two groups, young 1000 / 500. Block 1:
`copy_index` group 0 -> group 1 with factor 2.0 (=2000); block 2: `copy`
(all indices) with factor 0.5 -> group 0 = 500, group 1 = 1000. Final
values read by generic `target_item`s: exact unit test of both variants
and both factors.
