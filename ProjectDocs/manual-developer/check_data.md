# check_data

## Files and functions

- `database.cc:5021` — `check_data_integrity( void )`. Loops over all data items
  (`idat<MDAT`), and for each one with `data_required[idat]>0` that is
  `external[idat]` it iterates over the active indexes (`db_max_index` +
  `db_active_index`). If the item is active at index `index` but the required
  item `ireq = data_required[idat]` is not (`!db_active_index(ireq,index,...)`),
  it prints a `Warning` and counts it. A final summary line is printed with the
  total number of missing items.
- `initia.cc:95` — global `long int check_data=-NO;` (default: check off).
- `top.cc:129` — reads the keyword once at startup with
  `db( CHECK_DATA, 0, &check_data, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS )`
  and calls `check_data_integrity()` right after (`top.cc:135-136`) when
  `check_data==-YES`.
- `database.cc:255-258` — keyword registration:
  `strcpy(name[CHECK_DATA],"check_data")`, `type = INTEGER`, `data_length = 1`,
  `no_index[CHECK_DATA] = 1`.
- Enum `CHECK_DATA` in `tochnog.h:159` / `tochnog-mod.h:152` (must stay in sync).

## Implementation details

- The check only runs when `check_data==-YES`; the value is read through the
  `ival` argument because the type is `INTEGER`, not `dval`.
- "Required" means `data_required[idat]` holds the keyword of the companion
  data item that must be active at the same index.
- It relies on `db_name()` for the human-readable item names and on the
  existing `db_active_index` / `db_max_index` accessors.
- It does not abort the run; it only reports. Warnings are printed with
  `cout` directly (not through `pri()`).

## External dependencies

None. Uses only `db_*` accessors and globals already in scope in `database.cc`.

## Hardcoded parameters / pending refactorings

- Output goes to `std::cout`; routing the warnings through `pri()` would make
  them honor `check_warning -no`.
- It only covers items whose `data_required` is set. Strong validation of
  combinations is performed elsewhere by `check()` (e.g. `CHECK_COMBINATION`);
  `check_data` is complementary and only reports the missing companions.
- No per-item limit or early exit: it always scans all `MDAT` items.
