# control_print_history_factor

## Archivos y funciones

- `print_hi.cc` → `print_history()` (`print_hi.cc:22`) — the whole
  feature lives in the header of this function.
- `database.cc:881-886` — keyword registration:
  `strcpy(name[CONTROL_PRINT_HISTORY_FACTOR], "control_print_history_factor")`,
  type `DOUBLE_PRECISION`, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0`, class `CONTROL`, required for `CONTROL_PRINT_HISTORY`.
- Enum `CONTROL_PRINT_HISTORY_FACTOR` in `tochnog.h` / `tochnog-mod.h`.

## Detalles de implementación

- At `print_hi.cc:44-50` the code checks `db_active_index(
  CONTROL_PRINT_HISTORY_FACTOR, icontrol, VERSION_NORMAL )`. If active it
  allocates `factor_d`, initializes it to 1.0 with
  `array_set( factor_d, 1., DATA_ITEM_SIZE )`, then reads the user factors
  with `db( CONTROL_PRINT_HISTORY_FACTOR, icontrol, idum, factor_d, ldum,
  VERSION_NORMAL, GET )`.
- Per data item the printed value is multiplied by
  `( factor_d ? factor_d[idat] : 1. )` before writing to the `.his` file
  (same pattern as `print_data_versus_data`).
- `factor_d` defaults to 1.0 for all items when the keyword is absent, so
  existing inputs behave unchanged.

## Dependencias externas

- Core database accessors (`db()`, `db_active_index()`), `array_set()`
  from the utility layer. No external library.

## Parámetros hardcodeados / refactorizaciones pendientes

- `factor_d` is indexed by the *data item position* in the print record,
  which couples it to the order of `CONTROL_PRINT_HISTORY`. If the two
  lists drift out of sync the factors apply to the wrong items silently.
  A refactor could key factors by data item name/index instead of position.
- The "default to 1.0" initialization is duplicated in `print_hi.cc` and
  `print_da.cc`; a shared helper for "read factor list or default to 1.0"
  would remove the duplication.
- `DATA_ITEM_SIZE` caps the number of factors; inputs with more data items
  than this are truncated silently.
