# control_print_dof_line_time

## Implementación

- Stored as a CONTROL INTEGER record (`-yes`/`-no`). Read in
  `print_dof_line_point()` of `print_dl.cc`; `time_current` is fetched
  with `db( TIME_CURRENT, 0, ... )`.
- Written as the first line of each file of the call:
  ```
  if ( time==-YES )
    out << "# time " << time_current << "\n";
  ```
  (gnuplot comment format, manual 6.280).

## Diseño / decisiones

- The comment is written once PER CALL per file (the file is opened in
  append mode, so with several time steps the file contains one
  `# time` line per call, each followed by that call's data lines).
  "The first line of each file" (manual) is interpreted per written
  block; the alternative (only the very first line of the whole file)
  would make the time series useless.
- Format decision: `# time <value>` with `TN_PRECISION`.

## Detalles

- The comment is written even if no point of the call was found (the
  file then contains only the comment).

## Pendiente

- None.
