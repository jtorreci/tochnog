# condif_heat_edge_normal (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enums `CONDIF_HEAT_EDGE_NORMAL` and
  its 10 companions (headers in sync).
- `database.cc` — 11 registrations (master DOUBLE, 1 value; companions
  cross-required to the master, same shapes as the
  `groundflow_flux_edge_normal` family).
- `check.cc` — all 11 require the `condif_temperature` unknown.
- `area.cc` — the edge-normal machinery is SHARED with the groundflow
  water flux: `MTYPES 6 -> 7`, `type[6] = CONDIF_HEAT_EDGE_NORMAL`,
  `type_area[6] = CONDIF_HEAT_EDGE_NORMAL_GEOMETRY`. The three
  groundflow-specific blocks (element/element_group/element_side
  restrictions; sine/time load; the nodal application with factor and
  node restrictions) were generalized with two static helpers:
  `flux_edge_is_master()` and `flux_edge_companion(master, which)`
  mapping each master to its 9 companion enums. The application targets
  `temp_indx/nder` for the heat master (vs `pres_indx` for groundflow).
  The groundflow behaviour is bit-identical (regression:
  `groundflow_flux_edge` green).

## Verification

Test `condif_heat_edge` (suite 73/73): 2D column quad4 x2,
conductivity 0.1, q=0.1 injected at the bottom edge (geometry_line),
T=0 prescribed at the top: steady Fourier gives T(bottom)=2.0 and
T(middle)=1.0 — both targets exact.

## Pending

- `_element_side` uses the same pair convention as groundflow (skip if
  element not listed; the side number is informational for the filter).
- No dedicated test for `_node`/`_element_node`/`_element_node_factor`
  variants (same generic code path as groundflow, verified there).
