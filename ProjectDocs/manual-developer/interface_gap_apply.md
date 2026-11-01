# interface_gap_apply — developer

## Implementation

- **Enum**: `INTERFACE_GAP_APPLY` (`tochnog.h` / `tochnog-mod.h`), after
  `GROUP_INTERFACE_GAP`.
- **Record** (`database.cc`): INTEGER, length 1, `no_index 1`, CONTROL.
- **Consumer**: `interface_element()` (`interface.cc`): with `-no` the
  `GROUP_INTERFACE_GAP` read is skipped and the default gap 0 (C87)
  applies.

## Verification

C60/C87 geometry (four quad4 block pairs with a quad4 interface, pulled
0.01, slid 0.002, gap 1.e20): `-no` gives 0 / 0 on elements 3, 13, 23, 33,
equal to the Professional 02-08-2026 without the gap record; `-yes` gives
10 / +-2 = the Professional with the gap. The Professional rejects the
keyword itself. Tests: `interface_gap_apply1`, `interface_gap_apply2`.

## Pending

- `control_interface_gap_apply` (mentioned in manual 6.796, no entry).
