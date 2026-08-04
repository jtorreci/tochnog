# materi_plasti_maximum_iterations

## Archivos y funciones

- `stress.cc` → the block that sets `max_plasti_iter` inside the
  per-integration-point stress routine (`stress.cc:460-479`). After the
  built-in defaults are assigned, a user value overrides them.
- `database.cc:2502-2506` — keyword registration:
  `strcpy(name[GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS],
  "group_materi_plasti_maximum_iterations")`, type `INTEGER`,
  `data_length = 1`, class `MATERI`, required for `GROUP_TYPE`.
- Enum `GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS` in `tochnog.h` /
  `tochnog-mod.h`.

## Detalles de implementación

- Defaults before the override (`stress.cc:460-465`):
  - `materi_plasti_f_nonlocal && viscoplasti` → `max_plasti_iter = 2`,
  - `viscoplasti` → `max_plasti_iter = 10`,
  - otherwise → `max_plasti_iter = MAX_ITER`.
- The override block (`stress.cc:470-478`) reads with
  `db( GROUP_MATERI_PLASTI_MAXIMUM_ITERATIONS, gr, &max_plasti_iter_user,
  ddum, length_it, VERSION_NORMAL, GET_IF_EXISTS )`. If present and the
  value is `< 1` it prints an error and exits; otherwise
  `max_plasti_iter = max_plasti_iter_user`.
- The loop guard is at `stress.cc:538`:
  `while ( !plasti_found && total_plasti_iter<max_plasti_iter )`.

## Dependencias externas

- Core database `db()` accessor only; no external library.

## Parámetros hardcodeados / refactorizaciones pendientes

- IMPORTANT GOTCHA: `get_group_data()` does NOT return INTEGER values — it
  writes into `idum`. Do not use it here; the implementation correctly uses
  `db()` directly with the group index `gr`.
- The magic numbers `2` and `10` for the viscoplastic defaults are
  hardcoded; consider named constants alongside `MAX_ITER`.
- The "at least 1" validation (`pri` + `exit`) does not carry the group
  number in the message; adding `gr` would make the diagnostic actionable.
- `max_plasti_iter` and `total_plasti_iter` are state reset per element in
  `stress.cc:69`; a refactor could group them into a small struct to avoid
  the scattered reset.
