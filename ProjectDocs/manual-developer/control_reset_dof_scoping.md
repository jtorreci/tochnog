# control_reset_dof (per-index application + multi-dof records)

## Implementation

Both changes live in data.cc, in the `CONTROL_RESET_DOF` block of
`data()`.

1. **Per-index application**:
   ```c
   if ( ireset!=icontrol ) continue;
   ```
   right after `db_active_index( CONTROL_RESET_DOF, ireset, ... )`. The
   semantic was verified against the Professional binary: a reset
   sharing the timestep index is rejected ("Error detected for data
   item : control_reset_dof, record : 20"), so the reset must run at its
   own control index. Without this, re-applying e.g. the void ratio
   reset at every timestep desynchronized the dof from the model history.

2. **Multi-dof iteration**: the record length is captured from the GET
   (`reset_dof_length = ldum`) and the per-dof work (the
   `materi_displacement_relative` resync, the value-constant block, the
   value-dof block and the spatial distribution block) is wrapped in
   ```c
   for ( idof_list=0; idof_list<reset_dof_length; idof_list++ ) {
     idof_reset = reset_dof[idof_list];
     ...
   }
   ```
   The node selection filters (`control_reset_geometry`/`_node`/
   `_element_group`) are computed once per record, before the dof loop.

## Verification

- The Professional rejects a reset at the timestep index (apply-once
  semantics confirmed empirically).
- Corpus: 99 -> 112 PASS with the batch (the multi-dof fix unlocks
  multi-dof resets in many tests; hypo2/hypo4 pass).
- hypo1/hypo3 of the corpus stay RUNFAIL: with the CORRECTED reset
  semantics the wolfersdorff kernel (hypo.c port) deviates 0.019%
  (hypo1, sigyy -862.766 vs -862.92) and 36% (hypo3, anisotropic
  initial state) from the Professional targets - a kernel calibration
  issue, not a reset issue (documented in SEGUIMIENTO).

## Pending

- The multi-dof iteration changed the behaviour of the tests that
  previously relied on the first-dof-only bug (hypo1/hypo3). The kernel
  calibration needed to bring them back within tolerance is pending.
