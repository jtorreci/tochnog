# control_print_beam_force_moment_switch

## Implementación

- Registered in `database.cc`:
  `CONTROL_PRINT_BEAM_FORCE_MOMENT_SWITCH` (INTEGER, length 1,
  `data_required = CONTROL_PRINT_BEAM_FORCE_MOMENT`).
- Read inside `print_beam_force_moment()` (print_beam_force_moment.cc):
  `-YES` sets `factor = -1.`, `-NO` keeps `factor = +1.`, anything else
  is `db_error`. The factor multiplies all 12 components after they are
  built and before the snap-to-zero.

## Diseño / decisiones

- Default (record absent) is no inversion; `-no` is accepted explicitly
  although the manual only defines `-yes` (documented difference).
- The distance column is NOT inverted (the switch applies to the 12
  force/moment components only, per the manual text).

## Pendiente

- None.
