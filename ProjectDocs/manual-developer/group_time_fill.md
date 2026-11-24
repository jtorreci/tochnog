# group_time_fill — developer

## Files and functions

- `database.cc` — `GROUP_TIME_FILL` (DOUBLE_PRECISION, 3, class
  GROUP_TIME, requires GROUP_TYPE); enum after `GROUP_TIME` in `tochnog.h`
  / `tochnog-mod.h`.
- `group.cc`:
  - `group_time_fill_factor(element_group)`: 1 without the record,
    otherwise the clamped linear factor at `time_current + dtime`.
  - `group_time_alive()`: without `GROUP_TIME`, takes birth_empty / death
    of `GROUP_TIME_FILL`.
- `elem.cc` — `elem()`: the alive check also runs for a `GROUP_TIME_FILL`
  group.
- `bounda.cc` — `bounda_node_dormant()`: the dormant-node scan also runs
  when only `GROUP_TIME_FILL` records exist.
- `materi.cc`:
  - `get_materi_density()` returns the density times the factor (every
    caller: gravity, mass, `calcul_force`, `general.cc`).
  - the pore pressure added to the total stress (groundflow pressure unknown
    and static phreatic branches) is scaled through the pressure factor
    `gpf`.
- `plasti.cc` — `plasti_rule()` returns early, like
  `options_skip_plasticity -yes`, while the factor is below 1.

## Implementation details

Measured on the Professional 02-08-2026 (two stacked quad4, see the user
page): the density factor is evaluated at the end of the step; Mohr-Coulomb
is elastic until the step that reaches birth_filled; with a static phreatic
level at y = 2 the whole response of the filled element scales with f
(disy at y = 2: -0.00625 at f = 0.25; GNU without the pressure scaling
-0.0025).

The pressure coupling term of the element matrix (groundflow pressure
unknown) is not scaled, the same as `groundflow_pressure_factor`.

## History

- F5 T5 U7 (2026-11-19): added.

## Hardcoded parameters / pending refactorings

- After birth_filled the Mohr-Coulomb column of the probe differs from the
  Professional (-0.01531 vs -0.01526 at t = 3.5). This was already the case
  with a plain `group_time` (-0.01537 vs -0.01526), so it is not caused by
  this record.
- The groundflow pressure-unknown branch was measured on 2026-11-23
  (`group_time_fill_pressure1`: groundflow_pressure, a pressure gradient
  10 / 0, filling factors 0.25 / 0.75 / 1 at t = 1 / 2 / 3, with and
  without `groundflow_consolidation_apply -yes`): equal to the Professional
  02-08-2026 in disy, sigyy and pres. Without a phreatic level the total
  pressure is h - rho*g*z (materi.cc), so the prescribed pressure does not
  change the skeleton load in that probe.
