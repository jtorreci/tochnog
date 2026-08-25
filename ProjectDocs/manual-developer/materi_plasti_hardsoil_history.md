# materi_plasti_hardsoil_history

## Implementación

- **Initia**: new branch in the initia parser in `input.cc` (after
  `materi_plasti_cap1_history`). It sets the global flag
  `materi_plasti_hardsoil_history` (initia.cc) and allocates the dof:
  `sph_indx = unknown_indx; n = 1; dof_type = -MATERI_PLASTI_HARDSOIL_HISTORY`
  (basename `sph` in `database.cc`, SHARED with
  `materi_stress_pressure_history`).
- **Update**: `dof.cc` (`parallel_new_dof_diagonal`): the running
  maximum of |p| over the stress dofs is written into the sph dof,
  gated on `(materi_stress_pressure_history ||
  materi_plasti_hardsoil_history) && materi_stress` (the same block as
  the lot 6 sph update, condition extended). No equation in
  `general.cc` (the dof is driven only by this update, pattern of
  `materi_stress_pressure_history`).
- **Basename**: `sph` (same concept, same dof; the manual gives no
  basename).
- **check.cc**: `MATERI_PLASTI_HARDSOIL_HISTORY` requires
  `materi_stress` (the update reads the stress dofs).
- **Enum**: `MATERI_PLASTI_HARDSOIL_HISTORY` in `tochnog.h` /
  `tochnog-mod.h` (after `MATERI_PLASTI_CAP1_HISTORY`).

## Física

The manual (4.22): "The history variable abs(p) ... It contains the
maximum pressure history." Same concept as `materi_stress_pressure_history`
(4.50): both track the running maximum of |p| (`p = -sig_mean`,
positive in compression). The hardsoil elastic law uses it for the
E50/Eur (first loading vs unloading/reloading) switch.

## Gotchas

- The two initias share `sph_indx`; if both are present, the second
  parse overwrites the dof_type with its own value — harmless (the
  update and the basename are identical).
- The value read by the elastic block is `old_unknowns[sph_indx]` (the
  interpolated maximum at the START of the step); the dof is raised
  during the step by `dof.cc`, so the decision compares against the
  history EXCLUDING the current step (same reasoning as lot 6).
