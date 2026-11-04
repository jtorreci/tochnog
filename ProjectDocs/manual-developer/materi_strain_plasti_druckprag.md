# materi_strain_plasti_druckprag

## Implementación

- **Initia**: new branch in the initia parser in `input.cc` (after
  `materi_strain_plasti_diprisco`). Sets the global flag
  `materi_strain_plasti_druckprag` (initia.cc) and allocates the dof:
  `dpepp_indx = unknown_indx; n = 6;
  dof_type = -MATERI_STRAIN_PLASTI_DRUCKPRAG` (basenames
  `eppdrpxx eppdrpxy eppdrpxz eppdrpyy eppdrpyz eppdrpzz` in
  `database.cc`).
- **Integration**: the consolidated per-model plastic strain RHS block
  in `materi.cc` (same table mechanism as
  [materi_strain_plasti_cap](materi_strain_plasti_cap.md)).
- **general.cc**: inertia and convection entries for
  `-MATERI_STRAIN_PLASTI_DRUCKPRAG`.
- **check.cc**: `MATERI_STRAIN_PLASTI_DRUCKPRAG` self-check.
- **Enum**: `MATERI_STRAIN_PLASTI_DRUCKPRAG` in `tochnog.h` /
  `tochnog-mod.h`.

## Física

Manual (4.39): "The plastic strain eps_kl_plas specifically for the
drucker-prager model is added to the node_dof records." Same tensor as
`materi_strain_plasti` with a dedicated dof.

## Gotchas

- In the pure-shear validation the dedicated dof equals the generic epp
  EXACTLY (0.910768), confirming the consolidated RHS mechanism records
  the same strain; the value itself reflects the coupled stress-dof
  accumulation of the cutting-plane return (not the local integration
  point value).

## Model gating (C107, 2026-11-04)

The old consolidated block (h_lump weight, every group) was removed. The
dof is now one of the nine entries of the per-model block of
`materi.cc` (inside the strain recovery block, after the
total/elastic/plastic loop): `get_group_data(..., GET_IF_EXISTS)` on
`group_materi_plasti_druckprag` switches it on for the group, and the increment is
`volume * w_rec * inc_epp / dtime` with the `strain_recovery_weight` of
`materi_strain_plasti`. `general.cc` lists the type in the strain
recovery weight list, so the lumped denominator matches.
