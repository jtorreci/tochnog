# check_element_node

## Archivos y funciones

- `check.cc` → `check_element_node( long int check_element_node_switch )`
  (`check.cc:1068`) — iterates all elements and compares every pair of nodes;
  on a duplicate it prints
  `Error: element <ielem> has duplicate nodes.` to `std::cout` and calls
  `exit(TN_EXIT_STATUS)`.
- `input.cc:773-779` — called at the end of the data part, right after
  `end_data`. The switch defaults to `-YES` and is read with
  `db( CHECK_ELEMENT_NODE, 0, &check_element_node_switch, ddum, ldum,
  VERSION_NORMAL, GET_IF_EXISTS )` before calling the check.
- `database.cc:260-263` — keyword registration:
  `strcpy(name[CHECK_ELEMENT_NODE], "check_element_node")`, `type = INTEGER`,
  `data_length = 1`, `no_index[CHECK_ELEMENT_NODE] = 1`.
- Enum `CHECK_ELEMENT_NODE` in `tochnog.h:160` / `tochnog-mod.h:153`
  (must stay in sync).

## Detalles de implementación

- Early return: if the switch is `-NO`, or if the element table is empty
  (`db_max_index( ELEMENT, max_elem, VERSION_NORMAL, GET )` returns `max_elem < 0`),
  the function does nothing.
- For each active element it fetches `ELEMENT` and derives the node count as
  `nnol = length - 1` (first entry is the element type). It then runs a double
  loop `inol < jnol` over the nodes and aborts on the first `el[1+inol] ==
  el[1+jnol]` match.
- The switch is a LOCAL variable in `input.cc` defaulting to `-YES`, NOT a
  global: the default active state comes from the `-YES` initializer, not from
  `initia.cc`. If the keyword is absent, `GET_IF_EXISTS` leaves it at `-YES`.
- Allocation `get_new_int(1+MNOL)` is per call; freed with `delete[] el`.

## Dependencias externas

Uses the database layer (`db()`, `db_active_index()`, `db_max_index()`) and
the memory allocator `get_new_int()`. No external libraries.

## Parámetros hardcodeados / refactorizaciones pendientes

- The buffer size `1+MNOL` assumes a single element row never exceeds `MNOL`
  nodes; a `db` size check on `length` would make this defensive.
- Error output goes to `std::cout`; routing it through `pri()` / `tn.log`
  would be consistent with `exit_tn()` error reporting.
- `max_elem` and `nnol` are only used for iteration; a `while`-loop with
  `db_max_index` could be reused by a generic "validate all elements" helper.
- The switch is local to `input.cc`; other `check_*` features use a global
  read in `top.cc`, so the two patterns should be unified.
