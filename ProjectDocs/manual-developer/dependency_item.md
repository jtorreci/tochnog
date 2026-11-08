# dependency_item modifiers (developer)

## Where

- `tochnog.h` / `tochnog-mod.h`: `DEPENDENCY_GEOMETRY`, `DEPENDENCY_METHOD`,
  `DEPENDENCY_NUMBER`, `DEPENDENCY_TYPE` after `CONTROL_DEPENDENCY_APPLY`.
- `database.cc`: INTEGER records, class DEPENDENCY, `data_required`
  `DEPENDENCY_ITEM` (geometry 2 values, the others 1).
- `group.cc`, `get_group_data()`: the dependency loop.
  - Group match `dependency_item[1]==gr || ==-ALL`; with `-ALL` a group
    without its own record is skipped.
  - `dependency_geometry_contains()`: all element nodes inside the geometry
    (`geometry()`, `NODE_START_REFINED`, `PROJECT_EXACT`), the
    `element_geometry_present` rule; an invalid element index never matches.
  - The original record (DOUBLE records only) is read into `orig`; the int
    buffer of that `db()` call is `element_distribute` (full length, free
    at that point, audit contract).
  - Rows: `nrow = (length-n)/n`; one row broadcasts, otherwise one row per
    value (`nrow != norig` is a `dependency_diagram` error, measured on the
    Professional 02-08-2026). Without the original record the old behaviour
    stays (`nval` = fixed data length or `nrow`).
  - `dependency_interpolate()`: the old multi linear code, unchanged.
  - `dependency_combine()`: method/type (`atan`, `asin`, `acos` of the
    diagram value or of the value times `tan/sin/cos` of the original).
- `element_print_group_data_values_calculate( time_current )` (group.cc,
  called from `step_start`): real records now go through
  `get_group_data()` with the element node-average unknowns; the
  file-static `dependency_time_set` / `dependency_time` make the
  `-time_current` diagrams use the end time of the starting step (the time
  the element loop of that step uses). The element loop never sets them
  (serial print hook only).

## Measured on the Professional 02-08-2026

Probes in a two-element square (groups 0/1, mohr_coul in radians): -all,
multiply, number, type tangent/sinus/cosinus, one-row broadcast, row per
value, geometry_point radii 0.3 / 0.45 / 0.75 (only the element with all
nodes inside), errors for number beyond the record and a 2-row diagram on a
3-value record.

## Pending

- C122: the Professional evaluates `-time_current` at the middle of the step
  and integrates `-total_linear` incrementally; the GNU keeps the end of the
  step (analytic `E(t) eps`).
- `print_group_data` uses the element node average of the unknowns for dof
  dependencies (the Professional reports the last integration point).
- Space-coordinate dofs `-x/-y/-z` of `dependency_item` (manual 6.406) are
  not handled by the GNU dof lookup.
