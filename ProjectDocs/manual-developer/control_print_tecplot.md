# control_print_tecplot

## Implementation

- `top.cc`: `control_print_tecplot` other than `-no` calls
  `print_tecplot( ival )`; `print_tecplot_calculation` calls
  `print_tecplot( ival, 1 )` from `exit_tn` (`miscel.cc`).
- `print_te.cc` `print_tecplot()` picks the file first:
  `-SEPARATE_INDEX` -> `<data_file_base>_%04ld.plt` (ICONTROL),
  `-SEPARATE_SEQUENTIAL` -> `<data_file_base>_%04ld.plt` from the static
  `tecplot_sequential`; both are new files (`ios::trunc`). Otherwise
  `<data_file_base>.plt`: the static `tecplot_started` makes the first
  print of the run truncate it, later ones append. The end-of-run skip
  (`tecplot_printed` / `tecplot_time`) only tracks the `<base>.plt` prints.
- A new file forces `print_tecplot_new_mesh = -YES` (title, variables,
  coordinates and connectivity); appended zones share the connectivity of
  the first zone (`D=(1,2,FECONNECT)`).

## Verification (C99, 2026-11-02)

Professional 02-08-2026 probes: `-yes` -> `tey.plt` with every step, 5
states again after a second run in the same directory (the file is
emptied); `-separate_index` -> `tesi_0000.plt` (t 0.2 only),
`tesi_0012.plt` (t 0.5); `-separate_sequential` -> `tess_0000.plt` ..
`tess_0004.plt`. GNU tests `test-2026/control_print_tecplot_names1` (with
a stale `<base>.plt`), `_separate_index1`, `_separate_sequential1`.

## Hardcoded / pending

- 4-digit width for both separate switches (Professional).
- Zone layout: BACKLOG C27.
