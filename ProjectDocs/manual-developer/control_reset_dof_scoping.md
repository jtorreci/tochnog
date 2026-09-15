# control_reset_dof (per-index application + multi-dof records)

## Implementation

Both changes live in data.cc, in the `CONTROL_RESET_DOF` block of
`data()`.

1. **Per-index application**:
   ```c
   if ( ireset!=icontrol &&
        !db_active_index( CONTROL_TIMESTEP, ireset, VERSION_NORMAL ) )
     continue;
   ```
   right after `db_active_index( CONTROL_RESET_DOF, ireset, ... )`. The
   semantic was verified against the Professional binary: a reset
   sharing the timestep index is rejected ("Error detected for data
   item : control_reset_dof, record : 20"), so the reset must run at its
   own control index. Without this, re-applying e.g. the void ratio
   reset at every timestep desynchronized the dof from the model history.

   **Reconciliation (2026-09-15, delete3 regression)**: the extra
   `!db_active_index( CONTROL_TIMESTEP, ireset, ... )` clause keeps the
   legacy GNU every-step behavior ONLY for resets whose index also
   carries a timestep (a configuration the Professional rejects). For a
   reset at a non-timestep index the strict gate is preserved. See the
   section below for the two tests that pin the two branches.

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

## delete3 vs creset_geom reconciliation (2026-09-15)

Commit `f142bdf` removed the strict gate to fix `creset_geom`, which
re-opened the `delete3` regression: with the gate removed the reset at
index 5 fires at every timestep of the index-10 loop, re-injecting
`sigyy=1`, and the `post_point_dof 1 -sigyy` fixed point moves from
`0.600002285218` to `0.920015` (rc=1). Both tests are authoritative, so
the gate was refined instead of removed:

- **delete3** (`control_reset_dof 5 -sigyy` + `control_timestep 10`):
  index 5 has no timestep, so the reset is a ONE-SHOT initialisation and
  fires only at `icontrol==5`. Result restored to
  `0.600002285218` (target 0.6 +/-1e-4).
- **creset_geom** (`control_reset_dof 0` / `1` + `control_timestep 0` /
  `1`): both reset indices carry a timestep, so they keep the legacy
  every-step behavior. Without it, reset 1 (constant 0, geometry brick)
  does not fire at step 0; the brick `element_dof` IP history stays 0.5
  and step 1 solves `hisv0 = -1` instead of 0.
- **Diagnosis (instrumented firing trace)**: delete3 attempts the reset
  at `icontrol=5` (time 0) and at every `icontrol=10` step (0.2..1.0);
  creset_geom attempts both resets at `icontrol=0` (time 0.1) and
  `icontrol=1` (time 0.2). The gate must therefore depend on whether
  the reset's own index is a timestep index, not on `icontrol` alone.

The Professional rejects a reset sharing a timestep index, so the
legacy branch is unreachable for Professional-aligned inputs and the
sfnet-197 classification is unchanged (A/B byte-identical to HEAD,
187 PASS / 10 RUNFAIL / 0 PARSE; validation-suite 233 A/B identical
201/32/0; internal suite 16/16).

## Pending

- The multi-dof iteration changed the behaviour of the tests that
  previously relied on the first-dof-only bug (hypo1/hypo3). The kernel
  calibration needed to bring them back within tolerance is pending.
- A reset at a timestep index is still applied at EVERY step, including
  steps before its own index (legacy GNU). The Professional does not
  accept such inputs; if a future input needs "start at step i" scoping
  the condition must be revisited.
