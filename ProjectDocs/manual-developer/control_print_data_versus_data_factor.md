# control_print_data_versus_data_factor

## Archivos y funciones

- `print_da.cc` → `print_data_versus_data()` (`print_da.cc:22`) — the
  whole feature lives in the header of this function.
- `database.cc:814-819` — keyword registration:
  `strcpy(name[CONTROL_PRINT_DATA_VERSUS_DATA_FACTOR],
  "control_print_data_versus_data_factor")`, type `DOUBLE_PRECISION`,
  `data_length = DATA_ITEM_SIZE`, `fixed_length = 0`, class `CONTROL`,
  required for `CONTROL_PRINT_DATA_VERSUS_DATA`.
- Enum `CONTROL_PRINT_DATA_VERSUS_DATA_FACTOR` in `tochnog.h` /
  `tochnog-mod.h`.

## Detalles de implementación

- At `print_da.cc:42-48` the code checks `db_active_index(
  CONTROL_PRINT_DATA_VERSUS_DATA_FACTOR, icontrol, VERSION_NORMAL )`. If
  active it allocates `factor_d`, initializes to 1.0, and reads the user
  factors with `db( CONTROL_PRINT_DATA_VERSUS_DATA_FACTOR, icontrol, idum,
  factor_d, ldum, VERSION_NORMAL, GET )`.
- Per data item the value is multiplied by
  `( factor_d ? factor_d[idata] : 1. )` (`print_da.cc:70`) before appending
  to `tn.dvd`.
- Output file is opened as `ofstream out( "tn.dvd", ios::app )` with
  `out.precision(TN_PRECISION)` (`print_da.cc:53-55`).

## Dependencias externas

- Core database accessors, `array_set()` utility. No external library.

## Parámetros hardcodeados / refactorizaciones pendientes

- Identical structure to `control_print_history_factor`: `factor_d` is
  indexed by data item position in the `CONTROL_PRINT_DATA_VERSUS_DATA`
  record, so list ordering must stay in sync. Refactor to key factors by
  data item name/index, or extract a shared "read factor list or default
  to 1.0" helper used by both `print_da.cc` and `print_hi.cc`.
- The output filename `tn.dvd` is hardcoded at `print_da.cc:53` (consistent
  with upstream Tochnog, but could be a `control_print` option).
- `DATA_ITEM_SIZE` caps the number of factors silently.
