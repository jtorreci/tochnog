# control_groundflow_seepage_apply

## Where

- Enum: `CONTROL_GROUNDFLOW_SEEPAGE_APPLY` in `tochnog.h` +
  `tochnog-mod.h` (in sync).
- Registration: `database.cc` (class CONTROL, INTEGER x1).
- Combination check: `check.cc` — the three per-timestep groundflow
  control gates (`control_groundflow_consolidation_apply`,
  `control_groundflow_nonsaturated_apply`,
  `control_groundflow_seepage_apply`) require a groundflow pressure
  analysis.

## Bug fixed (combination check)

The three control records had `data_required = GROUNDFLOW`, which the
parser enforces as "an active GROUNDFLOW record at the SAME index as
the control". Group-type marker records are never stored, so EVERY
legitimate usage failed with "can only be used in combination with :
groundflow" (the marker only exists as a value inside `group_type`).
The data_required pairing was removed and the requirement is now
enforced in check.cc against the groundflow pressure DOF (like the
global no_index counterparts).

## Pending

- Consumption: the seepage branch of `bounda.cc` (seepage faces on
  `groundflow_seepage_geometry`/`_node`) does not read the control
  gate yet; tutorial_2 disables it for its linear stage.
