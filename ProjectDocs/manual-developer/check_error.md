# check_error

## Files and functions

- `pri.cc:22-30` — `pri( const char *s )`. The suppression lives at the top of
  the function: `if ( check_error==-NO && strstr( s, "Error" ) ) return;`
  (`pri.cc:25`), before the string is printed. Other `pri()` overloads are not
  affected.
- `initia.cc:92` — global `long int check_error=-YES;` (default: keep errors).
- `top.cc:126` — reads the keyword once at startup with
  `db( CHECK_ERROR, 0, &check_error, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS )`.
  The value arrives through the `ival` argument because the type is `INTEGER`.
- `database.cc:270-273` — keyword registration:
  `strcpy(name[CHECK_ERROR],"check_error")`, `type = INTEGER`,
  `data_length = 1`, `no_index[CHECK_ERROR] = 1`.
- Enum `CHECK_ERROR` in `tochnog.h:162` / `tochnog-mod.h:155` (must stay in sync).

## Implementation details

- The match is a substring test: any string containing "Error" is dropped when
  `check_error==-NO`, so messages such as `"Error: ..."` are filtered.
- Suppression happens before writing to `cout`, so the message is not flushed
  and does not reach `tn.log` output from `pri()`.
- Only messages that go through `pri( const char * )` are filterable.

## External dependencies

None. Uses only the global `check_error` and `strstr` from the C standard
library.

## Hardcoded parameters / pending refactorings

- LIMITATION: only messages routed through `pri( const char * )` are
  suppressed. Direct `cout` calls elsewhere in the code bypass the check.
- The filter is case-sensitive ("Error"): messages like `"error"` are not
  suppressed.
- Overloaded `pri()` variants (`const char*,const char*`, `int`, `long int`,
  `double`) do not apply the filter.
