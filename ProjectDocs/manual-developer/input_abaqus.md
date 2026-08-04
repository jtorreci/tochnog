# input_abaqus

## Archivos y funciones

- `input.cc` → `input_abaqus_read()` (`input.cc:1533`) — reads `abaqus.inp`
  line by line and writes `tochnog_abaqus.dat` with `node`, `element` and
  optional `group_materi_*` records.
- `input.cc` → global `input_abaqus_switch_global` (`input.cc:34`) — stores
  the `-yes/-no` switch of `input_abaqus`.
- `input.cc:535-553` — `idat==INPUT_ABAQUS` branch: parses the switch, saves
  it into `input_abaqus_switch_global`.
- `input.cc:556-576` — `idat==INPUT_ABAQUS_CONTINUE` branch: when the stored
  switch is `-YES` and the continue switch is `-YES`, calls
  `input_abaqus_read()`.
- `database.cc:2834-2870` — keyword registrations: `input_abaqus`,
  `input_abaqus_continue`, `input_abaqus_group`, `input_abaqus_mesh`,
  `input_abaqus_name`, `input_abaqus_set`.
- Enums `INPUT_ABAQUS`, `INPUT_ABAQUS_CONTINUE`, `INPUT_ABAQUS_GROUP`,
  `INPUT_ABAQUS_MESH`, `INPUT_ABAQUS_NAME`, `INPUT_ABAQUS_SET` in
  `tochnog.h` (`tochnog.h:625-630`) / `tochnog-mod.h` (must stay in sync).

## Detalles de implementación

- The reader parses `abaqus.inp` line by line with `in.getline(line,MCHAR)`,
  lowercases the line, trims leading spaces, then splits tokens with
  `strtok` by commas (`input.cc:1582-1583`).
- Section handling: `*node`, `*element` and `*elastic` keyword lines set the
  `in_nodes`, `in_elements`, `in_elastic` flags. The element type is taken
  from the `type=` attribute of the `*element` line (`input.cc:1594-1602`).
- The `*Elastic` material (Young's modulus and Poisson ratio) is captured and,
  when `input_abaqus_group==-YES`, written as
  `group_materi_elasti_young 0 <E>`, `group_materi_elasti_poisson 0 <nu>` and
  `group_materi_memory 0 -updated_without_rotation` (`input.cc:1711-1717`).
- `input_abaqus_name` and `input_abaqus_set` are read from the database
  (`db( INPUT_ABAQUS_NAME ... )`, `db( INPUT_ABAQUS_SET ... )`,
  `input.cc:1552-1567`) and applied as filters on each element
  (`input.cc:1682-1698`). `input_abaqus_group` is read with GET_IF_EXISTS
  (default `-YES`, `input.cc:1543, 1570-1571`).
- Type conversion chain (`input.cc:1629-1675`): compares the Abaqus type
  prefix and maps to a negated Tochnog element name and node count; any other
  type prints a warning and is skipped.
- Nodes are written with only `ndim` values (`input.cc:1619-1622`): the third
  coordinate is omitted when `ndim<3`. If all 3 values were always written,
  the extra value would break the subsequent parse.
- The generated file ends with `end_data` (`input.cc:1719`), so it can be
  included directly.

## Dependencias externas

None beyond the C++ standard library (`<fstream>`). Abaqus `.inp` parsing is
done manually; no Abaqus library is linked.

## Parámetros hardcodeados / refactorizaciones pendientes

- Input file `abaqus.inp` and output file `tochnog_abaqus.dat` are hardcoded
  (`input.cc:1573, 1579`).
- `input_abaqus_mesh` is registered but not implemented; the `-yes/-no`
  switch is never acted upon.
- Nset/Elset definitions are not converted: a missing `geometry_list` data
  record prevents writing group/geometry data. Set records are skipped
  silently.
- The type-conversion `else if` chain could be a table mapping Abaqus type
  prefixes to `{tochnog type, node count}`.
- `input_abaqus_switch_global` (file-scope global in `input.cc`) persists
  between the two data-part passes; it should be reset/owned by the
  `input_abaqus` branch to avoid stale state.
- Parsing relies on `strtok` (not thread-safe) and manual pointer trimming;
  a `std::istringstream`/`string_view` split would be cleaner.
