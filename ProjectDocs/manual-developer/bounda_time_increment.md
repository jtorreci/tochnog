# bounda_time_increment

## Files and functions

- `bounda.cc` — within `bounda()`: local `bounda_time_increment` (line 46),
  read of the keyword (lines 115–118, `GET_IF_EXISTS`) and use in the
  interpolation (lines 213–224).
- `tochnog.h` — enum `BOUNDA_TIME_INCREMENT` (line 144).
- `tochnog-mod.h` — mirror enum `BOUNDA_TIME_INCREMENT` (line 137).
- `database.cc` — keyword registration (lines 201–205).

## Implementation details

- Read: `db( BOUNDA_TIME_INCREMENT, iboun, idum, &bounda_time_increment, ... )`;
  only positive values activate the load-only mode
  (`if ( bounda_time_increment>0. )`).
- When active, `ninc = length_bounda_time` (NOT `length_bounda_time / 2`, which
  is the pair format branch). This keeps the interpolation loop intact while
  making every stored value a load.
- Interpolation (lines 213–224): `nload = length_bounda_time`; if
  `time_total >= bounda_time_offset`:
  `k = (long int) floor( (time_total - bounda_time_offset)
  / bounda_time_increment + 1.e-9 )`, clamped to `[0, nload-1]`, then
  `load = bounda_time[k]`. When `time_total < offset` the load stays 0
  (`found` remains 0).
- `database.cc`: `type = DOUBLE_PRECISION`, `data_length = 1`,
  `data_class = BOUNDA`, `data_required = BOUNDA_TIME`.

## External dependencies

None. Standard `floor` from the C library.

## Hardcoded parameters / pending refactorings

- The epsilon `1.e-9` added inside `floor` is a precision guard for `double`
  arithmetic (e.g. `0.03-0.02 == 0.0099999...`); the constant is hardcoded and
  not scaled to the increment magnitude — very small increments could still be
  off by one index.
- The branch duplication between the pair format and the increment format (both
  inside the `time` interpolation) could be unified.
- `bounda_time_increment` and `bounda_time_offset` share the same read pattern;
  the three modifier blocks (increment/offset, on_off, until_force) duplicate
  "read + validate + store" boilerplate that could be refactored.
