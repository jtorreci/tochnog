# materi_strain_plasti_camclay

## Implementation (F5 T2 V7)

- **Input:** `input.cc` parses the initia word, sets the flag
  `materi_strain_plasti_camclay` and the start index `camepp_indx`. It adds 6
  components of dof type `-MATERI_STRAIN_PLASTI_CAMCLAY`, with layout
  MATRIX.
- **Declarations:** the enum is appended in `tochnog.h` / `tochnog-mod.h`
  (F5 T2 V7 block), and the globals are in `initia.cc`.
- **Labels:** `database.cc` gives the labels `eppccxx..eppcczz` (dof label
  block) and the name.
- **Fill:** in `materi.cc`, inside the strain recovery block (the one that
  fills epe/epp/ept with `strain_recovery_weight`), after the block loop.
  The dof gets `volume * w_rec * inc_epp / dtime`, the same increment and
  weight as `materi_strain_plasti`. It is filled only when
  `get_group_data(GROUP_MATERI_PLASTI_CAMCLAY / _CAMCLAY_INCREMENTAL)` exists for the element group.
- **Recovery weight:** `general.cc` uses `strain_recovery_weight` for this
  dof type (recovery weight list), so the nodal value equals epp exactly.
  It also adds the type to the inertia list (inertia 1) and the
  convection list.
- **Other registrations:** `check.cc` (check_unknown) and `so.cc`
  (excluded from the linear shortcut).

## Verification

There is no Professional oracle (the Professional adds no dof). Test
`materi_strain_plasti_model3`: the camclay columns equal epp at every node (|diff| <= 1e-12), epp is
non-zero, and the camclay / mohr_coul / tension / vonmises columns of the
other models are 0 (`strain_plasti_model_check` in
`scripts/build_safe.sh`).

## Pending / notes

- The older per-model dofs (hardsoil, cap, compression, diprisco,
  druckprag, `ddacc43`) use a different mechanism. They are filled in
  every group whatever the model, with the h_lump weight, so on tria3
  their nodal values differ from epp. The Professional keeps them at 0
  under another model. Recorded as BACKLOG C107; it is not changed here
  because the change would alter existing output.
- Combined models are not split per yield surface.
