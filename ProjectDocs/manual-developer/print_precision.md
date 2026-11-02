# print_precision — developer

- Record `PRINT_PRECISION` (INTEGER, `no_index`). `TN_PRECISION`
  (`tochnog.h` / `tochnog-mod.h`) is now the global `tn_print_precision`
  (initia.cc, default 12), set in `top()` from the record (`db_error` below
  1). Every `out.precision(TN_PRECISION)` of the print routines follows it.
- `PROCESSORS_PARTITION` (`processors_partition`, 6.1006) is registered
  and not consumed: the GNU element loop gives each thread one static
  block of elements.
