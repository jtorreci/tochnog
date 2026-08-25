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
