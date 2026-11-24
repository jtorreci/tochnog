# control_input — developer

- Record `CONTROL_INPUT` (INTEGER, length 1, class CONTROL; `database.cc`
  before `control_materi_dynamic`, enum before
  `CONTROL_MATERI_DAMAGE_APPLY`).
- `input_control_input( icontrol )` (`input.cc`), called by `top()`
  (`top.cc`) at the entry of every control index, right after `ICONTROL`
  is set. Opens `<icontrol>.dat` in the working directory (missing:
  return).
- `input_simple_records( in, strict, mesh_read )` (`input.cc`, static):
  the simple record reader factored out of `input_runtime()` (no comments,
  defines, arithmetic or ranges). `strict` 0 keeps the runtime-file
  behaviour byte-identical (silent stop at the first bad word, `in >>
  index`, `any_runtime`); `strict` 1 (control_input) calls `check()` on
  each record, requires an integer index and stops with `Error for
  control_input detected.` / `Error detected near: <word>`.
- When a `node` or `element` record was read: zero `node_dof`,
  `node_start_refined` / `node_dof_start_refined` for the new nodes,
  zero `element_dof`, `element_dof_initialised` 0 and
  `nonlocal_element_info` for the new elements (the initialization of
  `top()`), then `mesh_has_changed( VERSION_NORMAL )`.
- Pending: `max_control` is computed before the loop, so a control record
  the file adds at a higher index than any other control record is not
  visited; the new-element integration-point state (BACKLOG C96).

## Differences with the Professional version

`CTRL-INPUT-NEWELEM` (class C), see the user page.
