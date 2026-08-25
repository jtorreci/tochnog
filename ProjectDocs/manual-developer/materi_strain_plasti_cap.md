# materi_strain_plasti_cap

## Implementación

- **Initia**: new branch in the initia parser in `input.cc` (after
  `materi_strain_plasti`). Sets the global flag
  `materi_strain_plasti_cap` (initia.cc) and allocates the dof:
  `capepp_indx = unknown_indx; n = 6;
  dof_type = -MATERI_STRAIN_PLASTI_CAP` (basenames
  `eppcapxx eppcapxy eppcapxz eppcappyy eppcapyz eppcapzz` in
  `database.cc`, same layout as `materi_strain_plasti`).
- **Integration**: `materi.cc` — one CONSOLIDATED RHS block drives ALL
  per-model plastic strain dofs (hardsoil + cap + compression +
  diprisco + druckprag): a small table of {flag, indx} pairs adds
  `element_rhside += volume*h*(inc_epp)/dtime` per component. The dof
  accumulates the SAME plastic strain increment `inc_epp` computed by
  the stress law driver (the per-model initias are registration
  aliases pointing to this single mechanism; `materi_strain_plasti`
  keeps its own separate block because of the temperature coupling).
- **general.cc**: inertia and convection entries for
  `-MATERI_STRAIN_PLASTI_CAP` (copies of the `-MATERI_STRAIN_PLASTI`
  entries, same as the other per-model strain dofs).
- **check.cc**: `MATERI_STRAIN_PLASTI_CAP` self-check (requires the
  initia string, like `MATERI_STRAIN_PLASTI`).
- **Enum**: `MATERI_STRAIN_PLASTI_CAP` in `tochnog.h`/`tochnog-mod.h`
  (before `MATERI_STRAIN_PLASTI_HARDSOIL`; both headers must stay in
  sync).

## Física

Manual (4.35): "The plastic strain eps_kl_plas specifically for cap
models is added to the node_dof records." Same tensor as
`materi_strain_plasti` (4.33) but with a dedicated dof so the cap
plastic strain can be post-processed independently.

## Gotchas

- The dof accumulates the plastic strain increment regardless of which
  law produced it (same behavior as `materi_strain_plasti`; the manual
  expects the cap model to be used alone, documented).
- 2D calculations still use the 6-component layout.
- The legacy `group_materi_plasti_cap` reads the GENERIC epp dof
  (`new_unknowns[epp_indx...]`); it still requires
  `materi_strain_plasti` (pre-existing GNU design, unchanged).
