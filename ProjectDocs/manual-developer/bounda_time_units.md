# bounda_time_units

## Files and functions

- `bounda.cc` — within `bounda()`: local array `bounda_time_units[2]`
  (line 48), read of the keyword and scaling of `bounda_time` in the
  `BOUNDA_TIME` read block (lines 116–126).
- `tochnog.h` — enum `BOUNDA_TIME_UNITS` (line 151).
- `tochnog-mod.h` — mirror enum `BOUNDA_TIME_UNITS` (line 144).
- `database.cc` — keyword registration (lines 241–245).

## Implementation details

- Read: `db( BOUNDA_TIME_UNITS, iboun, idum, bounda_time_units, ldum,
  ..., GET_IF_EXISTS )`, defaulting to `{1., 1.}`.
- Applied immediately after `BOUNDA_TIME` is read (and only in that branch,
  i.e. when `bounda_time` is given in data form). Iterates over all entries:
  even indices (`iu%2==0`) are times and are multiplied by
  `bounda_time_units[0]` (`factor_time`); odd indices are loads/lengths and
  are multiplied by `bounda_time_units[1]` (`factor_length`).
- Scaling is skipped when both factors are `1.` (the default).
- `database.cc`: `type = DOUBLE_PRECISION`, `data_length = 2`,
  `data_class = BOUNDA`, `data_required = BOUNDA_TIME`.

## External dependencies

None.

## Hardcoded parameters / pending refactorings

- The factors are only honored in the data-form branch of `BOUNDA_TIME`
  (line 112). The file-based branch (`BOUNDA_TIME_FILE`) and the user branch
  (`BOUNDA_TIME_USER`) ignore the keyword silently.
- No `data_length` check: if the user supplies a single value,
  `bounda_time_units[1]` stays at its initialized `1.` without warning.
- The odd/even index convention assumes the pair format `time load`; a
  `bounda_time` with a single value (length 1) never scales anything.
