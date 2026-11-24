# bounda_factor

## Files and functions

- `bounda.cc` — within `bounda()`: local array `bounda_factor[4]` (line 42),
  read of the keyword (lines 157–159, `GET_IF_EXISTS`, zeroed before with
  `array_set`) and application in the time branch of the node loop
  (lines 450–459 and 468).
- `tochnog.h` — enum `BOUNDA_FACTOR` (line 138).
- `tochnog-mod.h` — mirror enum `BOUNDA_FACTOR` (line 131).
- `database.cc` — keyword registration (lines 166–170).

## Implementation details

- Read: `db( BOUNDA_FACTOR, iboun, idum, bounda_factor, ... )` fills a fixed
  array `bounda_factor[4]`; coefficients not present are zeroed by
  `array_set( bounda_factor, 0., 4 )` before the read.
- Application (time branch, `else` of the `bounda_constant` test): if any
  coefficient is non-zero, the node coordinates are fetched with
  `db( NODE, inod, idum, coord_start, ... )` and
  `load_factor = a0 + a1*x + a2*y + a3*z`, using only the active dimensions
  (`ndim>=1`, `ndim>=2`, `ndim==3`). The new value is
  `new_node_dof[iuknwn] = factor * load * load_factor` (line 468).
- When all coefficients are zero the node-coordinate lookup is skipped
  (`load_factor = 1.`).
- `database.cc`: `type = DOUBLE_PRECISION`, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0`, `data_class = BOUNDA`.

## External dependencies

None. Internal database API only.

## Hardcoded parameters / pending refactorings

- The array is hardcoded to 4 coefficients for 3D (`a0..a3`); no dimension
  check is done on the input length (`data_length = DATA_ITEM_SIZE` accepts any
  number of values).
- `bounda_factor` and `bounda_factor_parabolic_x` are mutually exclusive at the
  value level: the parabolic block runs after the linear one and overwrites
  `load_factor` if any of its coefficients is non-zero. If both are given, the
  parabolic factor silently wins — should be documented or rejected.
- The `db( NODE, ... )` coordinate lookup runs per boundary node per boundary;
  for large models this can dominate the boundary loop — could be hoisted or
  cached.
