# print_database_calculation + print_gid_calculation

## Implementation

- Keywords registered in `database.cc`: `PRINT_DATABASE_CALCULATION`
  and `PRINT_GID_CALCULATION` (INTEGER, no_index, global switches).
  Enums appended in `tochnog.h`/`tochnog-mod.h` (same enum, same
  order).
- Consumption in `exit_tn()` (`miscel.cc`): before the final
  `print_database(-1, VERSION_NORMAL, -YES)` and `print_gid(-YES)`
  calls the switches are read (GET_IF_EXISTS); `-no` skips the
  corresponding dump. Default `-yes` keeps the previous behavior of
  every run.
- Verified: `large1.dat` parses and runs (the corpus still bounds it
  by the 45 s timeout; see SEGUIMIENTO for the solver-bound timing
  analysis).

## Pending

- The intermediate per-control outputs (`control_print_database`,
  `control_print_gid`) are not gated by the global switches (manual
  semantics: only the final "calculation" dump).
