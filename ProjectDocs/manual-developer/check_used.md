# check_used

## Archivos y funciones

- `database.cc` → `check_used_report()` (`database.cc:4998`) — scans all
  data items and prints the ones defined but never read.
- `database.cc` → `db()` / `db_dbl()` / `db_int()` — mark the read flag:
  `db_read[data_number] = 1` is set inside the GET / GET_IF_EXISTS code
  paths (`database.cc:4237, 4245, 4548, 4571, 4694, 4716`).
- `miscel.cc` → `exit_tn()` (`miscel.cc:472`) — calls
  `check_used_report()` at the end of the run when `check_used == -YES`.
- `initia.cc` (`initia.cc:90`) — global `long int check_used=-NO;`.
- `top.cc` (`top.cc:124`) — reads the keyword once at startup with
  `db( CHECK_USED, 0, &check_used, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS )`.
- `database.cc:312` — keyword registration:
  `strcpy(name[CHECK_USED], "check_used")`.
- Enum `CHECK_USED` in `tochnog.h` / `tochnog-mod.h` (must stay in sync).

## Detalles de implementación

- `check_used_report()` iterates `idat` from `0` to `MDAT-1` and reports an
  item when three conditions hold simultaneously:
  - `external[idat]` — the item was defined in the data file, and
  - `max_index[idat][VERSION_NORMAL] >= 0` — it has at least one index, and
  - `!db_read[idat]` — it was never read during the calculation.
- The report prints the keyword name and a final count to `std::cout`
  (not `pri()`), so it appears in the standard output at exit.
- `db_read[]` is a plain `MDAT`-sized global in `database.cc:41`. It is
  only ever set to `1`, never cleared after reading, which is correct for
  the "used at least once" semantics.

## Dependencias externas

None beyond the core database layer. Uses `db_name()`, `external[]`,
`max_index[][]` globals from the database module.

## Parámetros hardcodeados / refactorizaciones pendientes

- The output goes to `std::cout`; on runs that suppress the standard output
  the report may be invisible. Consider routing it through the existing
  log mechanism (`pri()` / output file) or a dedicated `.log` file.
- The three conditions in the loop are ad-hoc; extracting them into a named
  predicate would make the "defined but unused" rule testable in isolation.
- `check_used` is a global read once in `top.cc`; it could live in the
  database module next to `db_read[]` instead of `initia.cc`/`top.cc`.
