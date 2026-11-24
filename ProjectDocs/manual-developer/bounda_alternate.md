# bounda_alternate

## Files and functions

- `bounda.cc` — within `bounda()`, inside the iboun loop (lines 94-111).
  - Locals declared at lines 46-47: `bounda_alternate_list[DATA_ITEM_SIZE]`,
    `bounda_alternate_n`, `iteration`.
  - Current iteration read: `db( NUMBER_ITERATIONS, 0, &iteration, ddum,
    ldum, VERSION_NEW, GET_IF_EXISTS )` (line 82).
  - List read: `db( BOUNDA_ALTERNATE, 0, bounda_alternate_list, ddum,
    bounda_alternate_n, VERSION_NORMAL, GET_IF_EXISTS )` (lines 99-100);
    `bounda_alternate_n` receives the number of listed bounda indices.
  - Rotation logic (lines 101-111): if `iboun` is found in the list at
    position `jalt` and `iteration % bounda_alternate_n == jalt`, the record
    should be skipped for that iteration.
- `database.cc` — keyword registration (lines 160-164): name
  `bounda_alternate`, `type = INTEGER`, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0`, `data_class = BOUNDA`.
- `tochnog.h` — enum `BOUNDA_ALTERNATE` (line 138).
- `tochnog-mod.h` — mirror enum `BOUNDA_ALTERNATE` (line 131), must stay in sync.

## Implementation details

- The record is always stored/read at index 0. The first input value after the
  keyword is the record index; the remaining values are the bounda indices, so
  `bounda_alternate_list[0]` is the first listed bounda index and
  `bounda_alternate_n` is the list length.
- The iteration counter comes from `NUMBER_ITERATIONS` (`VERSION_NEW`), which
  `top.cc` writes at each iteration (line 330); the same record is used by
  `contact.cc` and `stress.cc`.
- Intended rotation: with n listed indices, at iteration `i` the index in list
  position `i % n` is omitted. E.g. `bounda_alternate 0 10 20 30` (n=3) omits
  10, 20, 30, 10, ... in successive iterations.
- `data_required` is unset for `BOUNDA_ALTERNATE`, so the record works
  standalone and does not force a `bounda_dof`/`bounda_unknown` to exist.

## External dependencies

None. Internal database API only; reuses the existing `NUMBER_ITERATIONS`
record.

## Hardcoded parameters / pending refactorings

- BUG: the `continue` at bounda.cc:108 is inside the inner
  `for ( jalt=0; jalt<bounda_alternate_n; jalt++ )` loop, so it only skips the
  rest of the inner loop and does NOT skip the bounda application. The bounda
  is still applied every iteration. To actually omit the boundary, the skip
  must act on the outer `iboun` loop (e.g. via a flag checked after the search,
  or by restructuring the loop).
- The list size is bounded by `DATA_ITEM_SIZE`.
- The record index is hardcoded to 0 in the `db()` call (line 99); a nonzero
  input index would be silently ignored.
- The `BOUNDA_ALTERNATE` read happens for every `iboun` inside the loop; it
  could be hoisted out of the iboun loop since the list does not change.
