# groundflow_apply — developer

- Record `GROUNDFLOW_APPLY` (INTEGER, length 1, class CONTROL, no index;
  `database.cc` before `control_materi_dynamic`, enum before
  `GROUNDFLOW_CONSOLIDATION_APPLY` in `tochnog.h` / `tochnog-mod.h`).
- `groundflow_apply_transform()` (`general.cc`), called once in `top()`
  (`top.cc`) after the node / element initialization and before the start
  version is stored, when the record is `-no`:
  - zeroes the groundflow dofs (`dof_type` `-GROUNDFLOW_PRESSURE`,
    `-GROUNDFLOW_VELOCITY`, `-GROUNDFLOW_SATURATION`,
    `-GROUNDFLOW_PRESSURE_GRADIENT`) of every `node_dof`;
  - `general_remove_groundflow_labels()`: removes their dof labels from
    `bounda_dof` / `bounda_unknown` (a record left without dof is deleted
    with its `bounda_time`) and `control_reset_dof` (deleted when empty);
  - `general_delete_prefix()`: deletes the `groundflow_*` (except
    `groundflow_apply`), `group_groundflow_*`,
    `group_interface_groundflow_*` and `control_groundflow_*` records;
  - sets the file-static flag read by `groundflow_apply_off()`.
- `element()` (`elem.cc`): the type loop skips `-groundflow`
  (`general()` and `groundflow()`) when `groundflow_apply_off()`.
- Pending: records added later by `control_input` or the runtime file are
  not rewritten (the switch is applied once, at the start).
