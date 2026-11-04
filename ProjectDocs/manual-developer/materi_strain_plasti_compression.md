# materi_strain_plasti_compression

## Implementación

- **Initia**: new branch in the initia parser in `input.cc` (after
  `materi_strain_plasti_cap`). Sets the global flag
  `materi_strain_plasti_compression` (initia.cc) and allocates the
  dof: `cepp_indx = unknown_indx; n = 6;
  dof_type = -MATERI_STRAIN_PLASTI_COMPRESSION` (basenames
  `eppcmpxx eppcmpxy eppcmpxz eppcmpyy eppcmpyz eppcmpzz` in
  `database.cc`).
- **Integration**: the consolidated per-model plastic strain RHS block
  in `materi.cc` (same table mechanism as
  [materi_strain_plasti_cap](materi_strain_plasti_cap.md); the dof
  accumulates the same `inc_epp` as `materi_strain_plasti`).
- **general.cc**: inertia and convection entries for
  `-MATERI_STRAIN_PLASTI_COMPRESSION`.
- **check.cc**: `MATERI_STRAIN_PLASTI_COMPRESSION` self-check.
- **Enum**: `MATERI_STRAIN_PLASTI_COMPRESSION` in `tochnog.h` /
  `tochnog-mod.h`.

## Física

Manual (4.36): "The plastic strain eps_kl_plas specifically for the
compression model is added to the node_dof records." Same tensor as
`materi_strain_plasti` with a dedicated dof.

## Gotchas

- Same as the other per-model strain initias: the dof records the
  plastic strain increment of whatever law yielded (documented; the
  manual expects the compression model alone).

## Model gating (C107, 2026-11-04)

The old consolidated block (h_lump weight, every group) was removed. The
dof is now one of the nine entries of the per-model block of
`materi.cc` (inside the strain recovery block, after the
total/elastic/plastic loop): `get_group_data(..., GET_IF_EXISTS)` on
`group_materi_plasti_compression` or `group_materi_plasti_compression_direct` switches it on for the group, and the increment is
`volume * w_rec * inc_epp / dtime` with the `strain_recovery_weight` of
`materi_strain_plasti`. `general.cc` lists the type in the strain
recovery weight list, so the lumped denominator matches.
