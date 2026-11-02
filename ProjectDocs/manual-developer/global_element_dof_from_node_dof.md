# global_element_dof_from_node_dof (developer)

- Enum `GLOBAL_ELEMENT_DOF_FROM_NODE_DOF` (`tochnog.h` / `tochnog-mod.h`,
  before `LAST_DUMMY`), registered in `database.cc` next to
  `OPTIONS_ELEMENT_DOF`: INTEGER, length 1, `no_index`.
- `top.cc` `top()`, the pre-time-loop element loop that creates the zero
  `ELEMENT_INTPNT_DOF` and `ELEMENT_DOF_INITIALISED` records: when the
  element had no `ELEMENT_INTPNT_DOF` record, `options_element_dof==-YES`
  and the switch is `-NO`, `ELEMENT_DOF_INITIALISED` is written as 1, so
  `elem.cc` restores the (zero) integration-point states instead of
  interpolating them from the nodes (the `element_dof_initialised==0`
  path). Otherwise the record is 0 as before.
- Default `-YES` hardcoded in `top.cc` (`element_dof_from_node_dof`): the
  historical GNU behaviour (`element_dof_initialised` 0 =
  interpolation from `node_dof`). Manual and Professional default `-no`:
  open owner decision ELEMENT-DOF-FROM-NODE-DOF-DEFAULT (flip the
  initial value to `-NO`; about 200 suite inputs hold `node_dof` records).
- Not touched: `CONTROL_DATA_INITELDOF_GEOMETRY` (`data.cc`) still resets
  `ELEMENT_DOF_INITIALISED` to 0 for a later re-initialisation from the
  nodes; new elements of the data part (`input.cc`), `generate.cc` and
  `create.cc` keep 0.

Verification: Professional 02-08-2026 probes f5/v5b `d_def` = `d_no` (no
initial stress taken), `d_yes` (taken; = GNU before the record),
`d_yes_ano` / `d_no_ano` / `d_def_ano` (apply -no: nodal values used,
switch without effect), `d_yes_aonly` (Professional aborts with -only).
Tests `global_element_dof_from_node_dof1/2`.
