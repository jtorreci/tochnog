# mpc_node_number / mpc_node_factor — developer notes

Manual Professional 6.874/6.875. Consumed in the new file `mpc.cc`
(`mpc_node_apply()`), called from `bounda.cc::bounda()` after the bounda
records (also when there are no bounda records: the call sits after the
`end_of_bounda:` label).

## Semantics (verified empirically against the Professional 25-10-2023)

- The slave dof is marked `NODE_BOUNDED` (known quantity) and its value
  is recomputed EVERY equilibrium iteration from the current master
  values (`NODE_DOF` VERSION_NEW, i.e. the values of the previous
  iteration, the step-start copy on iteration 1).
- The solver (`so.cc`) excludes bounded rows from the global system: the
  slave's equilibrium equation is dropped, the master's equation is
  untouched. NO static condensation and NO force redistribution —
  measured: a free master + rigid link makes the Professional's own
  pardiso fail on a singular system; an external force on a slave is
  kept in the slave's `node_rhside` but contributes nothing to the
  solution (GNU matches: the slave row is excluded from the solve).
- Cosmetic .dbs difference: the Professional zeroes the ELEMENT
  contributions in the slave's `node_rhside` (only the external force
  remains); the GNU assembles the full element contributions there. The
  solution values are identical (both exclude the slave row).

## Implementation

`mpc_node_apply()` in `mpc.cc`:

1. Calls `mpc_linear_quadratic_generate()` (mesh-change detection
   inside, see `mpc_linear_quadratic.md` developer notes).
2. Called from TWO places: `bounda.cc::bounda()` after the bounda
   records (pre-solve), and `top.cc` after
   `parallel_sys_routine(&parallel_new_dof_diagonal)` (post-solve). The
   post-solve call re-syncs the slave with the freshly solved masters —
   without it the slave lags one solve behind when the masters are free
   unknowns (measured: mpc3 slave 0.060 vs the tie value 0.129 with 2
   iterations; with the sync the slave always equals the tie).
3. Loops over all active `MPC_NODE_NUMBER` records:
   - Layout `[node_0 dof_0 node_1 dof_1a dof_1b ... node_2 ...]`: node
     numbers are positive, dof keywords negative. The slave dof is
     resolved via `array_member(dof_label, ...)` (same as bounda).
   - Masters are collected; `MPC_NODE_FACTOR` values map 1:1 in the
     flattened master-dof order, missing factors default to 1.
   - The slave node must exist (`db_active_index(NODE, ...)`); deleted
     nodes (mesh deletion) are skipped, deleted master nodes contribute 0.

## Gotchas

- The records were previously registered in `database.cc` without
  consumption ("I don't know what to do" era → target failures).
- `mpc_apply` (6.859) is NOT consumed yet: records are always active.
  The keyword is not even registered — registering it + reading it is a
  3-line follow-up.
- `mpc_geometry`/`mpc_element_group` automatic generation remains
  PENDING (different machinery: geometry- and element-group-based
  searches, `locate`-style).

## Pending

- `mpc_apply` consumption.
- `mpc_geometry` family (mpc2 passes by coincidence: rigid translation
  of a single quad4 tied by geometry).
- `mpc_element_group` family (mpc7 rc=0 since 2026-09-02; mpc6 stays RUNFAIL — a condif single-field model whose node-3 temp reads 0 vs 0.5, own diagnosis pending).
- `mpc_geometry_method`/`mpc_geometry_tolerance` keyword registration.
