# materi_displacement_relative

## Implementación

- **Keyword del initia** in `input.cc`: sets `materi_displacement_relative=1`,
  `dis_rel_indx = unknown_indx`, registers the dof as `-MATERI_DISPLACEMENT_RELATIVE`
  (VECTOR, `ndim`). Validates (input.cc): requires `materi_displacement`,
  `materi_velocity`, `materi_velocity_integrated`.
- **Enum**: `MATERI_DISPLACEMENT_RELATIVE` (+ internal data item
  `MATERI_DISPLACEMENT_RELATIVE_REF`) in `tochnog.h` / `tochnog-mod.h`.
  Labels `disrx/disry/disrz` in `database.cc` (dof_label block).
- **Integración**: in `dof.cc`, alongside `dis_indx`:
  `node_dof_new[dis_rel] = node_dof[dis_rel] + node_dof_new[vel]*dtime`.
- **Referencia (cambio de timestep)**: in `top.cc`, inside the
  `control_timestep` increment loop: if the new `dtime_initial` differs from
  the persisted `MATERI_DISPLACEMENT_RELATIVE_REF`, all `dis_rel` dofs are
  reset to 0 and `_REF` is updated.
- **Referencia (reset de desplazamiento)**: in `data.cc`, inside
  `control_reset_dof`: when the reset dof resolves to `dis_indx`, all
  `dis_rel` dofs are reset to 0.

## Validación

- `mat_rel`: quad4 uniaxial with `control_timestep 0 0.1 1. 0.5 1.` (two dt
  stages). `disy`=2.0 (total), `disry`=1.0 (reset when dt changed → only the
  second stage).
- `mat_rel_reset`: quad4 with `control_reset_dof -disx` + constant 0 each
  step. `disy`=1.0 (accumulated), `disry`=0.1 (reset every step by the
  displacement reset → only the last step).
