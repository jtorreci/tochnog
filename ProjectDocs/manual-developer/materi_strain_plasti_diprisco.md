# materi_strain_plasti_diprisco

## Implementación

- **Initia**: new branch in the initia parser in `input.cc` (after
  `materi_strain_plasti_compression`). Sets the global flag
  `materi_strain_plasti_diprisco` (initia.cc) and allocates the dof:
  `depp_indx = unknown_indx; n = 6;
  dof_type = -MATERI_STRAIN_PLASTI_DIPRISCO` (basenames
  `eppdipxx eppdipxy eppdipxz eppdipyy eppdipyz eppdipzz` in
  `database.cc`).
- **Integration**: the consolidated per-model plastic strain RHS block
  in `materi.cc` (same table mechanism as
  [materi_strain_plasti_cap](materi_strain_plasti_cap.md)).
- **general.cc**: inertia and convection entries for
  `-MATERI_STRAIN_PLASTI_DIPRISCO`.
- **check.cc**: `MATERI_STRAIN_PLASTI_DIPRISCO` self-check.
- **Enum**: `MATERI_STRAIN_PLASTI_DIPRISCO` in `tochnog.h` /
  `tochnog-mod.h`.

## Física

Manual (4.37): "The plastic strain eps_kl_plas specifically for the di
Prisco model is added to the node_dof records." Same tensor as
`materi_strain_plasti` with a dedicated dof. The di Prisco law itself
(plasti.cc diprisco block) works with the SHARED history dof (hisv),
see [materi_plasti_diprisco_history](materi_plasti_diprisco_history.md).

## Gotchas

- In the axisymmetric diprisc1 rig the shear components of the dof stay
  exactly 0 (axisymmetry); the radial and hoop components are equal.
- `materi_plasti_diprisco_density` (12 history variables) is NOT
  implemented: the interpolation law between the loose and dense
  parameter sets is not documented in the manual (it only references
  external papers), so an implementation could not be verified
  (see SEGUIMIENTO-CONVERGENCIA.md).
