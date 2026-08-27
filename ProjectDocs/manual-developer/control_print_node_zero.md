# control_print_node_zero

## Implementación

- CONTROL INTEGER record (1 value: `-yes`/`-no`), read with
  `GET_IF_EXISTS` in `print_node()` (`print_node.cc`), pre-set to
  `-YES` (the manual default). `data_required = CONTROL_PRINT_NODE` in
  `database.cc`.
- Applied in the write loop, PER FILE (per selected part): when
  `zero == -NO` and the printed value of that part is zero, the line is
  omitted.

## Diseño / decisiones

- **Exact zero comparison** (`value == 0.`): the manual says "results
  with value zero"; an exact comparison keeps the semantics crisp and
  is fully deterministic for prescribed Dirichlet fields. Documented
  difference.
- The filter is per part, so a node can disappear from one file (its
  zero part) while remaining in another (its non-zero part). With a
  sorted file the filtered line drops out of the sorted sequence.
- `-yes`/`-no` validation: anything else -> `db_error`.

## Detalles

- Only the printed VALUE is tested, not the angle or the coordinates.

## Pendiente

- None.
