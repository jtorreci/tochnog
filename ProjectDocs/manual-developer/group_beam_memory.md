# group_beam_memory (developer)

## Implementation

- `beam.cc`, `beam_2d()` and `beam_local()`: `memory` starts at
  `-TOTAL_LINEAR` (default since F5 T2 V4d, 2026-11-02, owner decision
  BEAM-MEMORY-DEFAULT; `-UPDATED` before), the group record is read with
  `GET_IF_EXISTS`, `-TOTAL_LINEAR` is mapped onto
  `-UPDATED_WITHOUT_ROTATION` and any other value than `-UPDATED` /
  `-UPDATED_WITHOUT_ROTATION` stops with `db_error`.
- `-UPDATED_WITHOUT_ROTATION`: the initial coordinates are used for the
  beam axes and length. `-UPDATED`: old coordinates = initial + old
  displacement, new = old + vel * dtime; the beam axes follow the chord.
- `database.cc`: `GROUP_BEAM_MEMORY` (INTEGER, length 1, class BEAM).

## Tests

- `test-2026/group_beam_memory_default1` (no record: 0.24, Professional
  probe `mem_total_linear` element 3).
- `test-2026/group_beam_memory_updated1` (`-updated`: 0.2399880, the
  former default; the Professional stops on this input).

## Pending refactorings

- `-updated_linear` (Professional option) is not accepted.
