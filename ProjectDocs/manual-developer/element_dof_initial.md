# element_dof_initial

## Implementation

- Keywords registered in `database.cc`:
  - `ELEMENT_DOF_INITIAL` (DOUBLE_PRECISION, variable length per
    element index, class ELEMENT) — the initial past dof values.
  - `ELEMENT_DOF_INITIAL_APPLIED` (INTEGER, version_all=1, not
    external) — internal marker that the initial field was applied.
- Consumption in `general()` (`general.cc`), the routine that
  integrates the lumped inertia of the transient equations
  (`(new - old)/dtime`): at the element's birth step the "old" value of
  the inertia increment is replaced by the `element_dof_initial` value
  of that dof position. The marker is checked at VERSION_NORMAL; it is
  written once the birth step converged by `step_close` (`top.cc`,
  task==YES pass over the elements that carry the record), so EVERY
  Newton iteration of the birth step sees the same initial field and
  the following steps integrate from the regular step-old dofs.
- The dof value for position `j` is `values[min(j, n-1)]` (last value
  reused when fewer values than dofs are given); the same value applies
  to all element nodes.
- Verified against the Professional binary 25-10-2023: `phase1.dat`
  gives `node_dof 2 = 0.7` EXACT at t=1 and the same
  element/node inertia distribution (measured 2-step probe: step 2
  continues from 0.7 -> 0.85, proving the initial field only affects
  the birth step).
- Blast radius: only elements that carry an `element_dof_initial`
  record are affected (previously the keyword did not even parse);
  all condif/conduction tests of the corpus keep their results.

## Pending

- The per-node layout of the manual ("values for the dofs of all
  nodes") is not implemented: values are per-dof and uniform over the
  element nodes (the corpus tests use the per-dof single-value form).
- `element_dof_initial_specific_number` / `_specific_value` (6.423/
  6.424) are not registered.
