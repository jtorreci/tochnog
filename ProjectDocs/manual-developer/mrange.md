# mrange and mstring — developer

- Names only (initialization keywords): enums
  `MRANGE_MAXIMUM_RANGE_LENGTH`, `MSTRING_MAXIMUM_NUMBER_OF_STRINGS`
  (`tochnog.h` / `tochnog-mod.h` after `MATERI_WORK`; the short enum names
  are taken by the `MRANGE` / `MSTRING` macros), names `mrange` /
  `mstring` in `database.cc` after `materi_work` (so `db_number()` and the
  keyword probe accept them).
- Parsing: `input()` initialization loop (`input.cc`), branches `mrange`
  and `mstring` through `input_initia_limit()` (integer > 0, at most the
  compiled `MRANGE` / `MSTRING-1`, error naming the constant). The value is
  stored in `initialization_values`, so the `.dbs` repeats the record.
- Limits: file-static `mrange_limit` / `mrange_given`, `mstring_limit` /
  `mstring_given` in `input.cc`.
  - `input_mrange_check()` (only when `mrange` is given): after
    `range_scan()` of an index range with its words including the two
    `-ra` markers (`range_length`, = Professional), after `range_expand()`
    with the expanded length,
    and when a value range `-ra ... -ra` closes in the value loop with its
    word count (`value_range_words`).
  - The define reader (`input_read_string()`) stops at `mstring_limit`
    strings; without `mstring` the limit is `MSTRING-1` and the message is
    the original one (byte-identical).
- Hardcoded: `MRANGE` 500000 (`tochnog.h`, range buffers in `input.cc`,
  `bounda.cc`, `data.cc`, `filter.cc`, check in `range.cc`), `MSTRING` 800
  (`input.cc`, static `define_strings[MDEFINE][MSTRING]`).
- Pending: the run-time expansion of value ranges (`bounda.cc`,
  `data.cc`, `filter.cc`) is checked against `MRANGE` only (BACKLOG C95).

## Differences with the Professional version

See the user page: `MSTRING-PRO` (class D), `MRANGE-COMPILED` (class C);
BACKLOG C95 (value-range expansion against `mrange`).
