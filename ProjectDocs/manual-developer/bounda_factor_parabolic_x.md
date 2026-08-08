# bounda_factor_parabolic_x

## Files and functions

- `bounda.cc` — within `bounda()`: local array `bounda_factor_px[3]` (line 42),
  read of the keyword (lines 160–162, `GET_IF_EXISTS`, zeroed before with
  `array_set`) and application in the time branch of the node loop
  (lines 460–467 and 468).
- `tochnog.h` — enum `BOUNDA_FACTOR_PARABOLIC_X` (line 139).
- `tochnog-mod.h` — mirror enum `BOUNDA_FACTOR_PARABOLIC_X` (line 132).
- `database.cc` — keyword registration (lines 172–176).

## Implementation details

- Read: `db( BOUNDA_FACTOR_PARABOLIC_X, iboun, idum, bounda_factor_px, ... )`
  fills a fixed array `bounda_factor_px[3]`; coefficients not present are
  zeroed by `array_set( bounda_factor_px, 0., 3 )` before the read.
- Application (time branch, `else` of the `bounda_constant` test): if any
  coefficient is non-zero, the node coordinates are fetched with
  `db( NODE, inod, idum, coord_start, ... )` and
  `load_factor = a0 + a1*x + a2*x^2` (lines 464–466). The new value is
  `new_node_dof[iuknwn] = factor * load * load_factor` (line 468).
- When all coefficients are zero the node-coordinate lookup is skipped
  (`load_factor = 1.`).
- `database.cc`: `type = DOUBLE_PRECISION`, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0`, `data_class = BOUNDA`.

## External dependencies

None. Internal database API only.

## Hardcoded parameters / pending refactorings

- The array is hardcoded to 3 coefficients (`a0..a2`); no length check on the
  input (`data_length = DATA_ITEM_SIZE` accepts any number of values).
- `bounda_factor` and `bounda_factor_parabolic_x` are mutually exclusive at the
  value level: this block runs after the linear one and overwrites
  `load_factor` if any coefficient is non-zero. If both are given, the parabolic
  factor silently wins — should be documented or rejected.
- The `db( NODE, ... )` coordinate lookup runs per boundary node per boundary;
  for large models this can dominate the boundary loop — could be hoisted or
  cached.
- Only `x` is supported; a general parabolic factor in `y`/`z` would duplicate
  this block and could be generalized.
