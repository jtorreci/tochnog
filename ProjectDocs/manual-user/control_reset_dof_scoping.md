# control_reset_dof (per-index application + multi-dof records)

## Description

Two refinements of `control_reset_dof` (manual Professional 6.350):

1. **Per-index application**: the reset is applied only in the control
   step of ITS OWN index (`icontrol == ireset`), like every control
   record. The Professional REJECTS a reset that shares the index of the
   timestep (`control_reset_dof 20 ...` together with
   `control_timestep 20 ...` gives "Error detected for data item :
   control_reset_dof, record : 20" - measured with the Professional
   binary), so the reset runs ONCE in its own (timeless) control step
   before the timesteps.

   **Reconciliation with legacy GNU inputs (2026-09-15)**: the
   one-shot rule is exact for resets at a NON-timestep index. A reset
   whose index ALSO carries a `control_timestep` (a configuration the
   Professional rejects) keeps the legacy GNU behavior and is applied
   at every `data()` call, i.e. at every step. This is what the
   reconstructed `creset_geom` needs (reset index 1 shares the timestep
   index 1 and must also fire at step 0 so the brick `element_dof` IP
   history follows the node). See
   `manual-developer/control_reset_dof_scoping.md` for the mechanism and
   the delete3/creset_geom trade-off.

2. **Multi-dof records**: a single record may list several dofs
   (`-sigxx -sigyy -sigzz`); EVERY listed dof receives the reset value.
   Previously only the first dof of the record was reset (a latent bug:
   the record was read but only `reset_dof[0]` was applied).

## Usage

```
control_reset_dof               1  -sigxx -sigyy -sigzz
control_reset_value_constant    1  -100
control_reset_dof               2  -hyhis0
control_reset_value_constant    2  1.03

control_timestep                20  1.e-3 0.3
```

The resets at index 1 and 2 run once each (their own control steps);
the timesteps run at index 20.

## Notes

- All value variants (`_value_constant`, `_value_dof`, the spatial
  distributions) apply to every listed dof.
- dam_building uses a 17-dof reset (`-velix -veliy -velx -vely -eppxx
  ... -kapsh`); every listed dof is now reset.
- hypo2/hypo4 of the corpus pass with the per-index + multi-dof
  behavior.
