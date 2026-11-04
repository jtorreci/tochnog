# materi_strain_plasti_hardsoil

## Implementación

- **Initia**: new branch in the initia parser in `input.cc` (after
  `materi_strain_plasti`). It sets the global flag
  `materi_strain_plasti_hardsoil` (initia.cc) and allocates the dof:
  `hsepp_indx = unknown_indx; n = 6;
  dof_type = -MATERI_STRAIN_PLASTI_HARDSOIL` (basenames
  `epphsxx epphsxy epphsxz epphsyy epphsyz epphszz` in `database.cc`,
  same layout as `materi_strain_plasti`).
- **Integration**: `materi.cc` RHS block (after the `materi_strain_plasti`
  block): `element_rhside += volume*h*(inc_epp)/dtime` per component
  (the dof accumulates the plastic strain increment; also filled for
  `options_element_dof`). No old-value read is needed (the increment
  drives the dof, same as `materi_strain_plasti`).
- **general.cc**: inertia and convection entries for
  `-MATERI_STRAIN_PLASTI_HARDSOIL` (copies of the
  `-MATERI_STRAIN_PLASTI` entries).
- **check.cc**: `MATERI_STRAIN_PLASTI_HARDSOIL` self-check (like
  `MATERI_STRAIN_PLASTI`).
- **Enum**: `MATERI_STRAIN_PLASTI_HARDSOIL` in `tochnog.h` /
  `tochnog-mod.h` (after `MATERI_STRAIN_PLASTI`).

## Física

Manual (4.40): "The plastic strain eps_kl_plas specifically for the
hardsoil model is added to the node_dof records." Same tensor as
`materi_strain_plasti` (4.33) but with a dedicated dof so the hardsoil
plastic strain can be post-processed independently.

## Gotchas

- The dof accumulates the plastic strain increment regardless of which
  law produced it (same behavior as `materi_strain_plasti`; the manual
  expects the hardsoil model to be used alone, documented).
- 2D calculations still use the 6-component layout (the out-of-plane
  component is filled too).

## Model gating (C107, 2026-11-04)

The old consolidated block (h_lump weight, every group) was removed. The
dof is now one of the nine entries of the per-model block of
`materi.cc` (inside the strain recovery block, after the
total/elastic/plastic loop): `get_group_data(..., GET_IF_EXISTS)` on
`group_materi_plasti_hardsoil` switches it on for the group, and the increment is
`volume * w_rec * inc_epp / dtime` with the `strain_recovery_weight` of
`materi_strain_plasti`. `general.cc` lists the type in the strain
recovery weight list, so the lumped denominator matches.
