# include

## Archivos y funciones

- `input.cc` → globals `include_file_stream` (`std::ifstream`, `input.cc:32`)
  and `include_reading` (`long int`, `input.cc:33`) — hold the stream and the
  active flag for the file being included.
- `input.cc` → `input_read_string()` (`input.cc:853`) — reads the next token
  from `include_file_stream` instead of `cin` while `include_reading` is set
  (`input.cc:866-872`). If the included file ends without `end_data` it errors
  out with "Unexpected end of include file detected."
- `input.cc` → `input()` (`input.cc:40`) — the main data-part loop processes
  the included file's records. The `INCLUDE` branch (`input.cc:620-642`) opens
  the file, sets `include_reading=1`, and reads the first token.
- `input.cc:514-523` — when `end_data` is found while `include_reading` is
  active, the included stream is closed, `include_reading` is reset to `0`,
  and the main stream (`cin`) is restored.
- `database.cc:2821` — keyword registration:
  `strcpy(name[INCLUDE],"include")`.
- Enum `INCLUDE` in `tochnog.h` (`tochnog.h:622`) / `tochnog-mod.h` (must stay
  in sync).

## Detalles de implementación

- Nesting is blocked: `include_depth` (`input.cc:53`) is incremented to `1`
  when an include starts, and the `INCLUDE` branch errors out if
  `include_depth>=1` ("include files cannot contain an include",
  `input.cc:622-626`).
- The included file is NOT parsed by a separate loop. It is fed through the
  same `while ( strcmp(str,"end_data") )` loop (`input.cc:512`) because
  `input_read_string()` transparently switches the source stream. This is the
  intended design — the file only needs `end_data` at its end.
- IMPORTANT: do NOT add a duplicate sub-loop to parse the included file. The
  original implementation did that and it was fragile, producing the
  "I do not know: 0" error because the sub-loop consumed tokens the main loop
  still expected. The main loop must keep processing the included file.
- File names are read as raw strings (no comment skipping inside the file), so
  comments in an included file are not supported.
- `input_skip_comment()` is still applied after each read, which is why
  comments on `end_data`/blank handling follows the main input rules.

## Dependencias externas

None beyond the C++ standard library (`<fstream>`).

## Parámetros hardcodeados / refactorizaciones pendientes

- `include_depth` is a local variable in `input()` that only supports one
  level; deeper nesting would require a stack of streams.
- The comment restriction comes from the raw string reading; supporting
  comments inside included files would need `input_skip_comment()` to run on
  tokens read from `include_file_stream` too.
- `include_file_stream` and `include_reading` are file-scope globals in
  `input.cc`; they could be wrapped in a small struct so the include state is
  explicit and resettable.
