# nonlocal / nonlocal_name

## Implementation

- **`nonlocal` alias** (database.cc, `db_number()`): the short
  Professional name routes to `OPTIONS_NONLOCAL` (DOUBLE_PRECISION,
  length 1, `no_index=1`):
  ```c
  if ( !strcmp( str, "nonlocal" ) )
    return OPTIONS_NONLOCAL;
  ```
  The GNU machinery (nonloc.cc) already reads `OPTIONS_NONLOCAL` as the
  averaging radius.

- **`nonlocal_name` record** (database.cc, `db_initialize()`): new enum
  `NONLOCAL_NAME` (tochnog.h / tochnog-mod.h in sync, after
  `OPTIONS_NONLOCAL_SOFTVAR`); INTEGER, length 1, `no_index=1`. The
  value is a negative name (`-group_materi_plasti_...`), parsed by the
  standard INTEGER-with-name path of input.cc.

## Verification

- slope_nonlocal_refine parses the records and starts the nonlocal
  computation. The run aborts in the legacy neighbour search
  (nonloc.cc): the `NONLOCAL_ITEM_SIZE` (160) buffer of `node_nonlocal`
  overflows for dense meshes - a pre-existing jan-2014 bug, documented
  in SEGUIMIENTO (PENDIENTE).

## Pending

- The per-model gate of `nonlocal_name` (the GNU applies the nonlocal
  contribution to every model).
- The nonlocal neighbour-search buffer bug (NONLOCAL_ITEM_SIZE).
