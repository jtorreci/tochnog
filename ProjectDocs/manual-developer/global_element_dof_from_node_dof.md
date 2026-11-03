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
- Default `-NO` hardcoded in `top.cc` (`element_dof_from_node_dof`), as
  in the manual and the Professional (owner decision C100, 2026-11-03).
  Before C100 the default was `-YES` (`element_dof_initialised` 0 =
  interpolation from `node_dof`); 20 suite inputs that give their initial
  state in `node_dof` now write `-yes` explicitly.
- `data.cc` `data()`, `control_reset_dof` dof-list loop: with the switch
  `-NO` (read again there, default `-NO`) and `options_element_dof==-YES`,
  after the node values of a reset record are set, every continuum
  element (`element_intpnt_supported()`, no interface group) with
  `ELEMENT_DOF_INITIALISED` 1, an `ELEMENT_INTPNT_DOF` record and all its
  nodes inside the reset node filter gets, at every integration point,
  the reset value evaluated at the point coordinates (h of `pol()` times
  the node coordinates; index `ipoint*nuknwn+dof`; the dof-dependent part
  uses the point's own dof values). The per-point value is the lambda
  `reset_value_at()`, shared with the node loop. Without it the reset
  never reached the integration points (they no longer copy the nodes);
  the Professional resets them pointwise (`control_reset_value_exponent1`
  1.99093 / 2.08768, `control_reset_value_logarithmic_first1`,
  `hypo_masin_history_slots1` hyhis0 1.03). With `-YES` nothing changes:
  the first-step copy carries the reset.
- Pending refactoring: the switch is read in two places (`top.cc` and
  `data.cc`); a global variable set once in `top()` would avoid that.
- Not touched: `CONTROL_DATA_INITELDOF_GEOMETRY` (`data.cc`) still resets
  `ELEMENT_DOF_INITIALISED` to 0 for a later re-initialisation from the
  nodes; new elements of the data part (`input.cc`), `generate.cc` and
  `create.cc` keep 0.

Verification: Professional 02-08-2026 probes f5/v5b `d_def` = `d_no` (no
initial stress taken), `d_yes` (taken; = GNU before the record),
`d_yes_ano` / `d_no_ano` / `d_def_ano` (apply -no: nodal values used,
switch without effect), `d_yes_aonly` (Professional aborts with -only).
Tests `global_element_dof_from_node_dof1/2/3` (3 = default) and
`hypo_masin_history_slots1` (reset under the default).
