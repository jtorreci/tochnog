# group_materi_stress_null / _direction / _direction_automatic (developer)

## Where

- `database.cc`: `GROUP_MATERI_STRESS_NULL` (INTEGER list, DATA_ITEM_SIZE,
  not fixed), `GROUP_MATERI_STRESS_NULL_DIRECTION` (DOUBLE, 9),
  `GROUP_MATERI_STRESS_NULL_DIRECTION_AUTOMATIC` (INTEGER, 1), all MATERI;
  the words `xy`, `xz`, `yz` (enums `XY`, `XZ`, `YZ`); `check.cc` requires
  `materi_stress`.
- `materi.cc` `materi_stress_null()`, called right after `set_stress()`:
  builds the local frame (rows = axes; explicit record normalized, zero
  axis an error; automatic: x = `ELEMENT_NORMAL`, y = the least aligned
  global axis made perpendicular, z = x cross y; else identity), and
  applies the projection
  `t := t - sum_{nulled ab} t_loc_ab (a_a b_b + a_b b_a)` to `new_sig` and
  to every strain column of `ddsdde` (via `stress_indx`).
- `area.cc`: `element_normal_needed()` and `element_normal_automatic()`
  include the automatic switch, so `element_normal_set()` builds the
  normals and reports a zero normal of such a group.

## Verification

Professional 02-08-2026, homogeneous hex8: global -y, -xy, -x -yz,
rotated and non-unit axes, von Mises: identical to all printed digits.
2D boundary layer with the automatic frame: identical to the converged
Professional (30 iterations).

## Pending

- `new_sig_eff` (damage) is not projected.
- `global_normal` is not registered (BACKLOG C14).
