# materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils_history

Initialisation (manual Professional 4.21). Adds the two history variables of the bonded-soil Cam
clay to the `node_dof` records:

- `gencamhis0`: p_s, the stress memory of the unbonded soil;
- `gencamhis1`: p_m, the bonding (tensile strength p_t = k p_m).

No number follows the keyword. Give both initial values yourself, for example with
`control_reset_dof`. The model is described in
[group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils](group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils.md).

The history variables of all models share one block of `node_dof` records, so this initialisation
cannot be combined with `materi_history_variables` or another `materi_plasti_*_history`:

```
Error: materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils_history cannot be combined
with another history initialisation (materi_history_variables or materi_plasti_*_history).
```

Tests: `materi_plasti_gencam_bonded_iso1` (dofs), `materi_plasti_gencam_bonded_error2` (combination).
