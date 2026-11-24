# bounda_time_offset

## Files and functions

- `bounda.cc` — within `bounda()`: local `bounda_time_offset` (line 46), read
  of the keyword (lines 155–156, `GET_IF_EXISTS`) and use in the interpolation
  (lines 216–218).
- `tochnog.h` — enum `BOUNDA_TIME_OFFSET` (line 145).
- `tochnog-mod.h` — mirror enum `BOUNDA_TIME_OFFSET` (line 138).
- `database.cc` — keyword registration (lines 207–211).

## Implementation details

- Read: `db( BOUNDA_TIME_OFFSET, iboun, idum, &bounda_time_offset, ... )`.
- Used only in the load-only interpolation of `bounda_time_increment`: the
  index is `k = floor( (time_total - bounda_time_offset) / increment + 1.e-9 )`.
  If `time_total < bounda_time_offset` no load is applied yet (`load = 0`).
- Note: it has no effect in the pair format of `bounda_time` (the `else`
  branch), where the times come directly from the data.
- `database.cc`: `type = DOUBLE_PRECISION`, `data_length = 1`,
  `data_class = BOUNDA`, `data_required = BOUNDA_TIME`.

## External dependencies

None. Standard `floor` from the C library.

## Hardcoded parameters / pending refactorings

- The epsilon `1.e-9` inside `floor` was added to avoid double-precision errors
  when subtracting the offset (e.g. `0.03-0.02 = 0.009999999...`). Hardcoded;
  see the same caveat in `bounda_time_increment`.
- Keyword interaction is implicit: `bounda_time_offset` only makes sense with
  `bounda_time_increment`, but `data_required = BOUNDA_TIME` does not enforce
  the increment keyword — misuse is silently ignored.
- Shared read pattern with `bounda_time_increment`; could be grouped into one
  read block.
