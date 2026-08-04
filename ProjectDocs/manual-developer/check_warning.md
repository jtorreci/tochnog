# check_warning

## Files and functions

- `pri.cc:22-30` — `pri( const char *s )`. The suppression lives at the top of
  the function: `if ( check_warning==-NO && strstr( s, "Warning" ) ) return;`
  (`pri.cc:26`), before the string is printed. Other `pri()` overloads are not
  affected.
- `initia.cc:93` — global `long int check_warning=-YES;` (default: keep
  warnings).
- `top.cc:127` — reads the keyword once at startup with
  `db( CHECK_WARNING, 0, &check_warning, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS )`.
  The value arrives through the `ival` argument because the type is `INTEGER`.
- `database.cc:307-310` — keyword registration:
  `strcpy(name[CHECK_WARNING],"check_warning")`, `type = INTEGER`,
  `data_length = 1`, `no_index[CHECK_WARNING] = 1`.
- Enum `CHECK_WARNING` in `tochnog.h:174` / `tochnog-mod.h:167` (must stay in
  sync).

## Implementation details

- The match is a substring test: any string containing "Warning" is dropped
  when `check_warning==-NO`, so messages such as `"Warning: ..."` are filtered.
- Suppression happens before writing to `cout`, so the message is not flushed
  and does not reach `tn.log` output from `pri()`.
- For a feature's warning to honor this switch, it must print through `pri()`
  (not direct `cout`). `check_element_shape` was corrected accordingly in
  `polynom.cc:453-458`: the message is formatted with `snprintf` into a
  `char warn_str[MCHAR]` buffer and emitted with `pri( warn_str )`.

## External dependencies

None. Uses only the global `check_warning` and `strstr` from the C standard
library.

## Hardcoded parameters / pending refactorings

- Only messages routed through `pri( const char * )` are suppressed; direct
  `cout` calls elsewhere in the code bypass the check.
- The filter is case-sensitive ("Warning"): messages like `"warning"` are not
  suppressed.
- Many older features still print warnings directly with `cout` (e.g.
  `check_data_integrity` in `database.cc`); migrating them to `pri()` would
  make them honor this switch.
