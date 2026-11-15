# groundflow_phreatic_level — developer

## Files and functions

- `database.cc` — record `GROUNDFLOW_PHREATICLEVEL` (DOUBLE, variable,
  no_index 1, class GROUNDFLOW) + `GROUNDFLOW_PHREATICLEVEL_N` (the 3D table
  layout). The single level is the `groundflow_phreaticlevel` GNU record
  (registered since sfnet 2014) under its Professional name alias
  (`db_number`).
- `groundfl.cc`:
  - `groundflow_phreatic_coord()` — evaluates the level at a node
    (constant or `table_xy`/`table_xyz`) and returns location/static/total:
    `static = rho*g*(level - coord_vertical)`, `total = pres_dof + static`,
    with the static/total clamped to `pressure_atmospheric` (default 0).
  - `groundflow_phreatic_apply()` — free-surface/dry-zone bounding (below).

## Implementation details

- Free-surface condition (added 2026-09-05): when a SINGLE
  `groundflow_phreatic_level` is active (and no
  `groundflow_phreatic_level_multiple`), `groundflow_phreatic_apply()`
  bounds the pres dof to 0 at every node with
  `coord_vertical >= level - EPS_COORD` (constant or table level per node).
  Rationale: above the level the static part of `groundflow_phreatic_coord`
  is positive and gets clamped to the atmospheric pressure, so
  `p_total = pres_dof` there — bounding `pres_dof = 0` enforces the dry
  condition `p_total = 0` measured on the Professional (its head on/above
  the phreatic line is `h = rho*g*level`, i.e. `p_dynamic = 0` in the GNU
  split). The bound gives the saturated flow domain its free-surface
  Dirichlet at the level row, making the pressure solve well-posed (without
  it the mesh above the level is only connected to the bottom Dirichlet and
  the coupled solve breaks down / oscillates).
- Ordering: `bounda()` zeroes `NODE_BOUNDED`, calls
  `groundflow_phreatic_apply()` (free surface first) and applies the
  bounda_dof records afterwards, so explicit `-pres`/`-topres` bounds win.
- `groundflow_phreatic_bounda` (manual 6.570, alias in `db_number()`,
  F5 T5 U1b) is read first in `groundflow_phreatic_apply()`: its default is
  `-YES` (record absent included, measured on the Professional 02-08-2026:
  absent and `-yes` give the same `.dbs`); `-NO` skips the free-surface
  block, so only bounda_dof bounds remain. Tests
  `groundflow_phreatic_bounda1` (`-yes`) and `_bounda2` (`-no`).
- The legacy `groundflow_phreaticlevel_bounda` METHODS 1/2 (explicit head
  prescription `pres_dof = rho*g*level`) are untouched and take precedence
  over the automatic free surface when the record is present (they bound the
  same nodes with the same bounded flag).
- Inclined single level (F5 T5 U3b, 2026-11-15):
  `groundflow_phreatic_level_reference(z_ref)` returns the height of the
  first table point (`level[0]` for one value, `level[1]` in 2D,
  `level[2]` in 3D). `groundflow_phreatic_coord()` computes
  `p_total = pres_dof + rho*g*(z_ref - z)` (static pressure still local),
  the free-surface block bounds `pres_dof = rho*g*(z_head - z_ref)` with
  `z_head` the level above the node, and the `bounda_dof -topres`
  conversion (bounda.cc, level_source 1) uses `z_ref` instead of the local
  level. The element Darcy terms use `grad(pres_dof)`, so the dof must be a
  head shifted by one constant: bounding 0 under an inclined level gave a
  hydrostatic field without flow. Horizontal levels: `z_ref = level`,
  same results as before. Measured on the Professional 02-08-2026 (2x3
  quad4, level `0 1 2 2.5`, `groundflow_phreatic_inclined1`): the total
  pressures match to 1e-6.
- `groundflow_phreatic_project` (manual 6.582, record
  `GROUNDFLOW_PHREATIC_PROJECT`, INTEGER, no index): with `-yes` and a 2D
  level line, `z_head` = `groundflow_phreatic_project_height()`, the
  height of the nearest point over all segments (parameter clamped to
  [0,1]). Measured on the Professional: node (0,3) under (0,1)-(2,2.5)
  gets head -1.72 (point (0.96,1.72)); node (2,3) the end point 2.5. A 3D
  table with `-yes` stops (not measured on the Professional). Words other
  than `-yes`/`-no` are `db_error`.
- The multiple-level family is separate (per-element-group ownership +
  `_static`), intentionally NOT touched by this branch.

## Verification

- ground3 (level at the top of the mesh, no flow): unchanged rc=0.
- ground13 (no level; static-height region): rc=0 with node values identical
  to the Professional (pres_dof uniform, to_pres/st_pres/dy_pres).
- ground14/15/16 (level 0.25 + -topres): the dry zone is bounded (pres_dof 0
  above the level, matching the Professional), the bottom bound converts
  -10 to pres_dof -7.5 and the saturated zone converges to the linear
  dynamic-pressure field; the piping/lifting targets (2.3333/2.0) are
  reproduced exactly once the coupled consolidation state reaches the
  drained steady field (~2000 s in the GNU; see SEGUIMIENTO for the
  remaining 1 s-window discrepancy of the u-p coupling).
