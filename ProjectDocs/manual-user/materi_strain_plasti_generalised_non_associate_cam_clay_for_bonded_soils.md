# materi_strain_plasti_generalised_non_associate_cam_clay_for_bonded_soils

Initialisation (manual Professional 4.38). Adds the plastic strain tensor of the bonded-soil Cam clay
to the `node_dof` records, dofs `eppgencamxx`, `eppgencamxy`, `eppgencamxz`, `eppgencamyy`,
`eppgencamyz` and `eppgencamzz`.

Like the other per-model plastic strains (see `materi_strain_plasti_camclay`), it is filled only in
groups with
[group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils](group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils.md),
where it equals `materi_strain_plasti`; in the other groups it stays 0.

Test: `materi_plasti_gencam_bonded_iso1` (eppgencamyy equals eppyy).
