# check_memory_usage

## Archivos y funciones

- `miscel.cc` → `exit_tn( long int print_database_type )` (`miscel.cc:374`) —
  in the SAME block as `check_memory`, at the end of the function just before
  `db_close()` (`miscel.cc:475-487`).
- `initia.cc:97-98` — globals `long int check_memory_usage=-NO;` and
  `double check_memory_usage_result=0.;`.
- `top.cc:133` — reads the keyword once at startup with
  `db( CHECK_MEMORY_USAGE, 0, &check_memory_usage, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS )`.
  Note the value arrives through the `ival` argument because the type is
  `INTEGER`.
- `database.cc:280-288` — keyword registration:
  `strcpy(name[CHECK_MEMORY_USAGE],"check_memory_usage")`,
  `strcpy(name[CHECK_MEMORY_USAGE_RESULT],"check_memory_usage_result")`,
  `type = INTEGER` / `DOUBLE_PRECISION`, `data_length = 1`,
  `no_index[...] = 1` for both.
- Enum `CHECK_MEMORY_USAGE` and `CHECK_MEMORY_USAGE_RESULT` in
  `tochnog.h:165-166` / `tochnog-mod.h:158-159` (must stay in sync).

## Detalles de implementación

- Enabled when `check_memory_usage==-YES`.
- Uses `getrusage( RUSAGE_SELF, &ru )` and stores the peak RSS converted to GB
  (`ru.ru_maxrss / 1048576.`) in the global `check_memory_usage_result`
  (`miscel.cc:479-481`).
- The result is persisted with
  `db( CHECK_MEMORY_USAGE_RESULT, 0, idum_m, &check_memory_usage_result, mem_len, VERSION_NORMAL, PUT )`
  (`miscel.cc:484-485`).
- NOTA / gotcha: `mem_len` is a local `long int mem_len=1;` and MUST be 1.
  Passing `ldum=0` (as in the surrounding GET calls) makes `db` fail on PUT —
  a zero length is not a valid PUT length.
- Prints `Peak memory usage: ... GB.` to `std::cout` and then calls
  `cout << flush;` because `exit(0)` follows immediately.
- The block is shared with `check_memory`; enabling either flag runs it.

## Dependencias externas

- POSIX `getrusage` / `struct rusage` — requires
  `#include <sys/resource.h>` at the top of `miscel.cc` (`miscel.cc:23`).
- `db(..., PUT)` for persistence.

## Parameters hardcodeados / refactorizaciones pendientes

- The GB conversion hardcodes `1048576.` (KB→GB). On macOS `ru_maxrss` is in
  bytes, so the constant is wrong there; a platform-aware conversion
  (`#ifdef __APPLE__`) is pending.
- `mem_len=1` is duplicated as a magic value inline; it should be the keyword's
  `data_length` from `database.cc` (or a named constant) to prevent the PUT
  failure recurring.
- The peak-capture block is duplicated between `check_memory` and
  `check_memory_usage`; a shared helper that captures and optionally persists
  the peak would remove the duplication.
- Output goes to `std::cout`; routing it through `pri()` / `tn.log` would be
  consistent with the rest of Tochnog reporting.
