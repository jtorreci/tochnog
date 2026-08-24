# control_data_activate (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enum `CONTROL_DATA_ACTIVATE`.
- `database.cc` — registration: `INTEGER`, `DATA_ITEM_SIZE`,
  `fixed_length = 0` (variable list of item names + switch).
- `data.cc` — block in `data()`: reads the record, the switch is the
  LAST slot; for every listed item name with switch `-no`, all active
  indices are deleted via `db_delete_index` (deletion is the GNU
  "de-activation": `db_active_index` reports the records as gone and
  every consumer loop skips them). `-yes` is a no-op (input records are
  active by default; the deletion is destructive, so re-activation
  would need `control_data_put` — documented GNU difference).

## Consumer safety

Deleting an item mid-run is safe for consumers that check
`db_active_index` per record and derive their loop bounds from
`db_max_index` of the same item — the standard pattern in this codebase
(e.g. `bounda()` skips deleted `bounda_force` records). Orphaned
companion records (e.g. a `bounda_time` whose `bounda_force` was
deleted) are skipped the same way.

## Verification

Test `cda_activate` (suite 67/67): column pulled by `bounda_force`
(total +10 in vely) during block 0 -> disy(top) ~ +0.01 (sigyy = 10,
epsyy = 0.01, young 1000, poisson 0). Block 1 fires
`control_data_activate 1 -bounda_force -no`: the load disappears and the
quasi-static (`options_inertia -no`) equilibrium returns disy -> 0
(target 0.0 ± 4e-3). With the force still active the target would fail
(~0.01).

Model gotchas found while calibrating: `materi_displacement` requires
`derivatives` in initia AND `group_materi_memory -total_linear`
(`-updated_without_rotation` forbids displacement unknowns, stress.cc).
