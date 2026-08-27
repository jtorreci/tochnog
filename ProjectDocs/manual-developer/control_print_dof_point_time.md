# control_print_dof_point_time

## Implementación

- Stored as a CONTROL INTEGER record (`-yes`/`-no`). Read in the
  `is_line=0` branch of `print_dof_line_point()` in `print_dl.cc`
  (validated `-yes`/`-no`).
- Same behaviour as `control_print_dof_line_time`: `# time <time_current>`
  as the first line of each call's block (gnuplot comment format, manual
  6.283).

## Diseño / decisiones

- Identical to the line variant; the point family has no `_n`, so each
  call writes exactly one data line per dof file after the comment.

## Detalles

- The comment is written even when the point is not found (file with
  only the comment).

## Pendiente

- None.
