# bounda_sine

## Files and functions

- `bounda.cc` — within `bounda()`:
  - record read / `nfreq` / `time_start` / `time_end`: lines 296-307.
  - start/end gating (`found`): lines 573-578.
  - first consumer (prescribed `new_node_dof`, displacement / velocity /
    acceleration boundary): lines 799-810.
  - second consumer (`node_rhside`, force boundary): lines 1174-1184.
- `database.cc` — keyword registration `BOUNDA_SINE` (lines 280-284):
  `type = DOUBLE_PRECISION`, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0`, `data_class = BOUNDA`.
- `tochnog.h` / `tochnog-mod.h` — enum `BOUNDA_SINE` (mirrored, same order).

## Implementation details

- Record layout (Professional manual 6.34):
  `[start_time, end_time, freq_0, amp_0, freq_1, amp_1, ...]`.
- Read: `time_start = bounda_sine[0]; time_end = bounda_sine[1];`
  `nfreq = ( db_len( BOUNDA_SINE, iboun, VERSION_NORMAL ) - 2 ) / 2;`
- Evaluation (both consumer sites):
  `frequency = bounda_sine[2 + 2*i];`
  `amplitude = bounda_sine[3 + 2*i];`
  `value = amplitude * sin( 2*pi*frequency*time_total )`.
- The derivative (needed by the prescribed-displacement consumer) adds
  `amplitude * 2*pi*frequency * cos( 2*pi*frequency*time_total )`.
- Gating: `found = time_current >= time_start && time_current <= time_end;`
  so the dof is prescribed only inside the window and free outside it.
- `time_total = time_current + dtime` is the time used to evaluate the sine
  (end of step), while the gating uses `time_current` (start of step),
  consistent with the pre-existing start-time gate.

## Format correction

The GNU consumer previously used `nfreq = ( len - 1 ) / 2` and
`bounda_sine[1 + 2*i]` / `bounda_sine[2 + 2*i]`, i.e. it assumed
`[start_time, freq, amp, ...]` and ignored `end_time`. For the Professional
input `bounda_sine 20 0.0 1.e8 1.e1 1.e-3` it therefore read `freq = 1e8`,
`amp = 10` instead of `freq = 10`, `amp = 1e-3`, so the model ran effectively
unexcited. The fix aligns both consumer sites and the `nfreq` count with the
Professional layout.

`end_time` handling is implemented (not PENDING): manual 6.34 states "The sine
loads will be only imposed after start_time, and will not be imposed anymore
after end_time", and the `bounda_time` section states that a dof not
prescribed in such way is free and determined by the governing differential
equations. The gate therefore deactivates the prescription for
`time_current > time_end`.

## Pending

- Frequency-0 special case (manual 6.34: "setting a frequency to 0 enforces
  tochnog to use a constant static value of the specified amplitude") is not
  implemented. Currently a zero frequency contributes 0.
- No length-based disambiguation of the legacy GNU 3-value form
  (`start freq amp`) is provided. The blast-radius audit found no GNU input
  using that form, so it was not added. A legacy input would now be read as
  `nfreq = (3-2)/2 = 0` (zero load) instead of raising an error.

## External dependencies

None. `sin` / `cos` from the C library.

## Hardcoded parameters / pending refactorings

- `PIRAD` is the global pi constant.
- The record is read with `data_length = DATA_ITEM_SIZE` and
  `fixed_length = 0` (variable length), so the number of `(freq, amp)` pairs
  is inferred from `db_len`; there is no upper-bound check.
- `time_end` is a local of `bounda()`, reset per boundary inside the
  `BOUNDA_SINE` branch; it is not persisted in the database.
