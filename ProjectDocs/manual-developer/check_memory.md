# check_memory

## Archivos y funciones

- `miscel.cc` → `exit_tn( long int print_database_type )` (`miscel.cc:374`) —
  the block sits at the end of the function, after `check_used_report()` and
  just before `db_close()` (`miscel.cc:475-487`).
- `initia.cc:96` — global `long int check_memory=-NO;` (default: off).
- `top.cc:132` — reads the keyword once at startup with
  `db( CHECK_MEMORY, 0, &check_memory, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS )`.
  Note the value arrives through the `ival` argument because the type is
  `INTEGER`.
- `database.cc:275-278` — keyword registration:
  `strcpy(name[CHECK_MEMORY],"check_memory")`, `type = INTEGER`,
  `data_length = 1`, `no_index[CHECK_MEMORY] = 1`.
- Enum `CHECK_MEMORY` in `tochnog.h:164` / `tochnog-mod.h:157`
  (must stay in sync).

## Detalles de implementación

- Enabled when `check_memory==-YES`.
- Uses `getrusage( RUSAGE_SELF, &ru )` and reports `ru.ru_maxrss`
  (peak resident set size, in KB on Linux) converted to GB via `/ 1048576.`
  (`miscel.cc:479-481`).
- Prints `Peak memory usage: ... GB.` to `std::cout` and then calls
  `cout << flush;`.
- IMPORTANT: `exit(0)` is called right after this block (`miscel.cc:494`).
  Without the explicit flush the buffered output is lost, so the flush is
  required, not cosmetic.
- The block is shared with `check_memory_usage`; enabling either flag runs it
  (`if ( check_memory==-YES || check_memory_usage==-YES )`).

## Dependencias externas

- POSIX `getrusage` / `struct rusage` — requires
  `#include <sys/resource.h>` at the top of `miscel.cc` (`miscel.cc:23`).
  `RUSAGE_SELF` counts the calling process only (not children).

## Parameters hardcodeados / refactorizaciones pendientes

- The GB conversion hardcodes `1048576.` (KB→GB). On macOS `ru_maxrss` is in
  bytes, so the constant is wrong there; a platform-aware conversion
  (`#ifdef __APPLE__`) is pending.
- Output goes to `std::cout`; routing it through `pri()` / `tn.log` would be
  consistent with the rest of Tochnog reporting.
- The block is duplicated between `check_memory` and `check_memory_usage`;
  the peak-value capture could be factored into one helper called once.
