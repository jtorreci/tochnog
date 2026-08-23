# bounda_water

## Files and functions

- `bounda.cc` — inside `bounda()`, in the time branch of the per-dof loop.
  - Read per iboun: `db( BOUNDA_WATER, iboun, &bounda_water, ... )` with
    `GET_IF_EXISTS` (lines 203-205); `bounda_water = 0` by default.
  - Applied when `bounda_water==-YES && iuknwn==pres_indx` (line 518).
  - Direct computation `sp = fg[ndim-1] * dens * ( wl - coord[ndim-1] )`
    (line 541), using `force_gravity_calculate( fg )`, the first value of
    `groundflow_density` read with `db_dbl` (it is `no_index`), and the first
    value of `groundflow_phreaticlevel` (lines 521-540).
  - Result assigned with `new_node_dof[iuknwn] = factor * sp` (line 542),
    where `factor` carries the `bounda_time` multiplier.
- `database.cc` — keyword registration (lines 264-268): name `bounda_water`,
  `type = INTEGER`, `data_length = 1`, `data_class = BOUNDA`,
  `data_required = BOUNDA_UNKNOWN`.
- `tochnog.h` — enum `BOUNDA_WATER` (line 155).
- `tochnog-mod.h` — mirror enum `BOUNDA_WATER` (line 148), must stay in sync.

## Implementation details

- The branch sits in the `else` of the `bounda_constant` check, so
  `bounda_constant` takes precedence over the water pressure.
- The vertical gravity component is `fg[ndim-1]`; the vertical node coordinate
  is `coord[ndim-1]`, so the code is dimension-agnostic (1D/2D/3D).
- `bounda_water` is read with `GET_IF_EXISTS`, so it is inactive when no record
  exists for `iboun`.
- The pore pressure dof is the `pres_indx` dof; the record prescribes `-pres`.

## External dependencies

- `force_gravity_calculate()` (`force.cc`) — provides the gravity vector,
  including the optional `force_gravity_time` factor.
- `groundflow_density` (`database.cc:2005-2009`) — `no_index`, read with
  `db_dbl` because it has no index.
- `groundflow_phreaticlevel` (`database.cc:2011-2016`) — only the first value
  is used here; the full interpolation (`table_xy`/`table_xyz`) is not applied.
- `pres_indx`, `bounda_time` (`factor`), `db`/`db_dbl`/`db_active_index`.

## Hardcoded parameters / pending refactorings

- Only the FIRST value of `groundflow_phreaticlevel` is used; the table-based
  (2D/3D) phreatic interpolation that `groundflow_phreatic_coord()` supports is
  ignored here.
- `groundflow_phreatic_coord()` (`groundfl.cc:110`) computes the same static
  pressure but then CLAMPS `static_pressure` to `pressure_atmospheric`
  (default 0), which annuls the pressure. `bounda_water` deliberately bypasses
  that function and does the direct calculation, so the clamping is not applied
  (also `total_pressure` is untouched there for consistency).
- `groundflow_density` is only read when `db_active_index( GROUNDFLOW_DENSITY )`
  is true; otherwise `dens` stays 0 and the pressure becomes 0. No fallback
  density is used.
- Keyword is an INTEGER flag (`-yes`/`-no`) stored per iboun; could be merged
  with `bounda_unknown` handling instead of being a separate record.

## Known limitation: the `pressure_atmospheric` clamp in `groundflow_phreatic_coord()`

`groundfl.cc:323-325` applies an upper clamp:

```c
if ( static_pressure>=pressure_atmospheric ) static_pressure = pressure_atmospheric;
if ( total_pressure>=pressure_atmospheric ) total_pressure = pressure_atmospheric;
```

Semantics (see `groundflow_pressure_atmospheric.md` for the full write-up):

- This is a **cap/clamp**, not a gauge conversion: the code never subtracts
  the atmospheric pressure; it truncates any value at the
  `groundflow_pressure_atmospheric` threshold (`database.cc:2476`,
  `no_index=1`, default **0** if unspecified).
- With `force_gravity (0,-1)` the static pressure `fg[dens]*(wl-y)` is
  NEGATIVE below the phreatic level (compression) and POSITIVE above it
  (suction). The clamp therefore caps the SUCTION side: with the default 0,
  suction above the phreatic surface is annulled ("no suction" behaviour),
  while compression below the level passes through untouched. A positive
  threshold keeps suction up to that value.
- `bounda_water` bypasses this clamp entirely and uses the direct formula,
  so it returns the full hydrostatic pressure on both sides of the level.
  This is intentional but means `bounda_water` and
  `groundflow_phreatic_coord()` give different numbers above the phreatic
  level unless `groundflow_pressure_atmospheric` is configured
  consistently.

`groundflow_pressure_atmospheric` is documented and verified in its own
manual pages (tests `groundflow_pressure_atm` / `_def`). Remaining open
question: whether the clamp should also apply to `bounda_water` for
consistency with the unsaturated-soil model.
