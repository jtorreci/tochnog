# control_print_history_smooth

## Implementación

- **Output**: `print_history_smooth()` in `print_hi.cc` (same file as
  `print_history`). Invoked from the control loop in `top.cc`, right after
  `print_history` for the same `icontrol`:
  ```
  print_history( ival, nval );
  if ( db_active_index( CONTROL_PRINT_HISTORY_SMOOTH, icontrol, VERSION_NORMAL ) )
    print_history_smooth( ival, nval );
  ```
  It receives the same `ival[]`/`nval` as `control_print_history` (the
  data item records).
- **Keyword**: `control_print_history_smooth` (data_class CONTROL,
  data_length DATA_ITEM_SIZE, type INTEGER) registered in `database.cc`,
  with `data_required = CONTROL_PRINT_HISTORY`.
- **New enum**: `CONTROL_PRINT_HISTORY_SMOOTH` in `tochnog.h` /
  `tochnog-mod.h` (kept in sync).
- **Window sizes**: read per control record; `nsmooth` = number of integer
  values. If `nsmooth >= nset` (one value per history set), `smooth_d[iset]`
  is used per set; otherwise `smooth_d[0]` applies to all sets.
- **Moving average**: a persistent ring buffer per history set
  (`static double* smooth_buf[DATA_ITEM_SIZE]` etc.). The buffer is
  reallocated when the window size changes. Each call stores the current
  value, then averages the last `min(count, N)` values.

## Diseño / decisiones

- The buffer is `static` (persists across calls within one run). The set
  index `iset` maps 1:1 to the history sets of `control_print_history`, so
  the buffer is indexed by `iset`.
- The data value resolution (data_item_name/index/number, incl. the
  `number /= nder` for matrix dofs) is copied from `print_history` to keep
  both files consistent.
- File names follow `print_history` exactly, prefixed with `smooth`.

## Detalles

- The value can be INTEGER or DOUBLE_PRECISION depending on the data item;
  both are handled (cast to double for the average).
- Ring buffer index: `smooth_count[iset] % smooth`; the count is the total
  number of stored values (used to average fewer than N at the start).

## Pendiente

- The ring buffers are never freed (static arrays; fine for a single run,
  but not re-entrant across repeated initialization of the control loop).
- `db_error( CONTROL_PRINT_HISTORY_SMOOTH, ... )` used for out-of-range
  `number` (mirrors print_history).
