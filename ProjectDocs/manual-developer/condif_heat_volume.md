# condif_heat_volume (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enums `CONDIF_HEAT_VOLUME` and its 8
  companions; `condif()` signature EXTENDED with `coord_ip[]` (needed
  for the spatial factor at the integration point); declaration of
  `user_condif_heat_volume()`.
- `database.cc` — 9 registrations (master DOUBLE, 1; `_geometry`
  INTEGER 2; `_user` INTEGER 1; lists/factors/tables variable length,
  cross-required).
- `check.cc` — all 9 require `condif_temperature`.
- `condif.cc` — new block after absorption: loops over active
  `CONDIF_HEAT_VOLUME` records; restrictions by `_element`,
  `_element_group` (array_member) and `_geometry` (all element nodes in
  the geometry via `geometry()`, NODE_START_REFINED); heat from the
  record or from `user_condif_heat_volume` when `_user -yes`; load from
  `_sine` (start_time + freq/amp pairs, area.cc pattern) or `_time`
  (`force_time`); spatial `factor` via `force_factor` at `coord_ip`;
  contribution `element_rhside[i] += volume * h[i] * heat * load *
  factor`.
- `elem.cc` — call site updated with `coord_ip` (per-integration-point
  coordinates, already computed in the ipoint loop).
- `user.cc` — `user_condif_heat_volume` stub following the
  `user_viscosity` pattern: explicit error until the user programs it.

## Gotchas

- `group_condif_flow` in 2D takes TWO values (flow vector components);
  one value is an input-parse error (found calibrating the edge test).
- The steady-state trick of the tests: `group_condif_capacity 0 0.0`
  drops the storage term, so one timestep reaches the steady solution.

## Verification (suite 73/73)

- `condif_heat_vol`: 1D bar2 x2, k=1, T=0 both ends, S=1 restricted by
  `_element` to element 2: T(center) = 0.25 analytic (0.5 with both
  elements heated — the A/B proves the restriction).
- `condif_heat_vol2`: `_factor 0. 1.` gives S(x)=x: T(center) = 0.5
  analytic (polynomial factor path).
