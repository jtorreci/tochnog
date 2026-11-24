# ranges

## Where

- `range.cc` `range_scan()`: reads an index range after the `-ra` read
  by `input()` (`input.cc`, index part of the data records) into the
  compact form `-RA items... -RA`, where a block is
  `-FROM a -TO b [-STEP s]`. Static helpers `range_scan_error()`,
  `range_read_word()` (= `input_read_string()` + `input_skip_comment()`
  + lower case) and `range_read_integer()`.
- `range.cc` `range_expand()`: expands the compact form into the list of
  numbers. Used for index ranges (`input.cc`, store loop) and for value
  ranges at the point of use: `bounda.cc` (bounda_dof/force node lists),
  `data.cc` (`data_delete`, `data_put`, `data_list`, `data_arith`),
  `filter.cc` (`print_filter`).
- Value ranges are read by the generic record reader in `input.cc`
  (`-ra`, `-from`, `-to`, `-step` become `-RA`, `-FROM`, `-TO`,
  `-STEP` through `db_number()`); `in_value_range` toggles on each
  `-RA` value and makes a non-integer number an input error.

## Implementation

- Index ranges are fully validated in `range_scan()`: `from > to`,
  `step <= 0`, a missing `-to`, a non-integer word and an empty range
  stop the run; after a block (and its optional `-step`) the loop reads
  the next word, so blocks and list items can follow in any order.
- `range_expand()` validates the blocks of value ranges, which never
  pass through `range_scan()`: `-FROM` must be followed by a number,
  `-TO` and a number (a negative entry is a keyword), `from <= to`,
  `step > 0`. It no longer assumes `-TO` is at `j+2` without checking.
  The loop is a `while`, so an empty value range gives
  `range_length` 0.
- Error texts: `Error in range, a -from is larger than a -to` is the
  Professional text.
- `range_test_scan()` / `range_test_expand()` keep the `MRANGE` bound.

## Bugs fixed (F4 T6f, 2026-10-12)

1. `range_scan()` stored `STEP` instead of `-STEP`; `range_expand()`
   tests `-STEP`, so an index `-step` range expanded `from..to` and then
   the step and its value as indices (duplicate-record error).
2. `from >= to` rejected `-from a -to a`.
3. A `-from` after `-to n` (second block) was an error.
4. A descending value block expanded to nothing, silently.
5. `9.0` in a value range was read with `atoi`.
6. `cin >> ival` / `cin >> str` raw reads (unbounded word, no defines,
   include files or comments): BACKLOG C56 residual, now through
   `input_read_string()` and its bounded `input_word()`.

Verification: audit probes (Engram 2758, 32 forms x index/value) give
the same expansion as the Professional 02-08-2026 for every accepted
form and an error for every rejected one; tests
`validation-suite/test-2026/range_*.dat` (5 rc 0 with targets equal to
the Professional `.dbs`, 3 expected errors) in `scripts/build_safe.sh`.

## Hardcoded parameters

- `MRANGE` (tochnog.h): maximum compact length and expanded length.

## Pending refactorings

- Value ranges are validated only when expanded (`range_expand()` at
  the point of use), not at input time; an invalid value range in a
  record that is never used goes unnoticed.
- `range_scan()` takes `d`/`d_is_set` by value only to pass them to
  `input_read_string()`; arithmetic results are ignored there (as the
  Professional, which rejects arithmetic words in a range).

## Differences with the Professional version

- Error message wording and exit code (GNU 1, Professional 0).
- Integer values with a decimal point outside ranges are accepted by the
  GNU (`atoi`), rejected by the Professional (DIVERGENCES `INT-DOT`).
