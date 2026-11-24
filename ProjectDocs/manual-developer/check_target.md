# check_target

## Archivos y funciones

- `miscel.cc` → `exit_tn()` — evaluates the `target_item` / `target_value`
  records. When a target is not met and `check_target == -YES` it writes an
  error to `tn.log` and calls `exit(TN_EXIT_STATUS)` (`miscel.cc:437-455`);
  when `check_target == -NO` it writes a "Note" to `tn.log` instead and
  continues (`miscel.cc:457-465`).
- `initia.cc` (`initia.cc:91`) — global `long int check_target=-YES;`.
- `top.cc` (`top.cc:125`) — reads the keyword once at startup with
  `db( CHECK_TARGET, 0, &check_target, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS )`.
- `database.cc:302-305` — keyword registration:
  `strcpy(name[CHECK_TARGET], "check_target")`, `type = INTEGER`,
  `data_length = 1`, `no_index[CHECK_TARGET] = 1`.
- Enum `CHECK_TARGET` in `tochnog.h:170` / `tochnog-mod.h:163` (must stay in sync).

## Detalles de implementación

- The target check lives inside `exit_tn()` (`miscel.cc`), iterating the
  `target_item` / `target_value` records and comparing each against the
  referenced data item with a tolerance.
- The only difference in behaviour is the branch inside the `!correct` block:
  `-YES` writes an error to `tn.log` and exits; `-NO` writes a "Note"
  (`check_target -no`) to `tn.log` and keeps going. The message text is
  duplicated between both branches.
- The global is read once with `GET_IF_EXISTS`, so an absent keyword keeps the
  `initia.cc` default (`-YES`).
- `TARGET_ITEM` and `TARGET_VALUE` are registered in `database.cc:3839` and
  `database.cc:3845`.

## Dependencias externas

None beyond the core database layer (`db()`, `db_dbl()`, `db_int()`,
`db_name()`, `db_type()`) and the `TN_EXIT_STATUS` exit code.

## Parámetros hardcodeados / refactorizaciones pendientes

- The log message is written with `ofstream out( "tn.log", ios::app )` in
  `exit_tn()`, duplicating the file-open logic that other logging paths use.
  A shared helper for writing to `tn.log` would avoid the copy-paste.
- The "error" and "note" branches duplicate almost identical text formatting;
  they could share a single write with a different prefix/verbosity level.
- `check_target` is a global read once in `top.cc`; the value is never
  re-read, so a change mid-run has no effect (probably intentional).
