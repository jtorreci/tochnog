# mesh_activate_gravity_time

## Implementación

- **Función**: `mesh_activate_gravity_factor()` in `mesh.cc`. Returns the
  gravity activation factor (0 before activation, 1 after, interpolated) for
  an element. Called from `materi.cc` after `force_gravity_calculate()`, the
  gravity vector is multiplied by the factor.
- **Keywords** (data_class CONTROL, index 0) in `database.cc`:
  `mesh_activate_gravity_time` (DOUBLE 2), `_element`/`_element_group`
  (INTEGER variable), `_geometry` (INTEGER 2), `_method` (INTEGER 1),
  `_stiffness_factor` (DOUBLE 1), `_time_initial` (DOUBLE 1),
  `_time_strain_settlement` (INTEGER 1), and `control_mesh_activate_gravity_apply`
  (INTEGER variable). Enums in `tochnog.h`/`tochnog-mod.h`.
- **Selección de elementos**: `_element` (range match), `_element_group`
  (group match), `_geometry` (all nodes inside). Without a selection record,
  all elements apply.
- **Interpolación**: the element activation interval is `time_start`..`time_end`
  (the manual interpolates by the lowest/highest element coordinate; here the
  full window is used for a single element — bottom-to-top interpolation is a
  simplification).
- **`_time_initial`**: before `time_of_birth` the factor is 0.

## Pendiente

- Nada pendiente del `mesh_activate_gravity_*`: method 1 y 2 implementados.
  Nota: con stiffness_factor bajo y dofs libres el modelo puede ser
  numéricamente inestable (documentado en el test mesh_act_grav2, que usa
  todos los nodos fijos).

## F4 T8 (2026-10-10, BACKLOG C46d)

- `mesh.cc` `mesh_activate_gravity_factor()` rewritten with three static
  helpers: `mesh_activate_gravity_selected(index, ...)` (selection
  records of one index; none = all), `mesh_activate_gravity_vertical()`
  (y in 2D, z in 3D, from `NODE_START_REFINED`), and
  `mesh_activate_gravity_range(index, zmin, zmax)` (vertical range of the
  selected elements, cached per index and per (time_current, dtime,
  max element)). The first index whose selection contains the element
  wins. Element window: `t_start_el = ts + (te-ts)(zmin_el-zmin)/(zmax-zmin)`,
  `t_end_el` with `zmax_el`. Time `time_current + dtime` (end of step;
  the start of step was worse against the Professional).
- Default method `-METHOD2` (manual 6.810; before `-METHOD1`). Before the
  element start time: factor 0, stiffness factor `stiff_factor_val`
  (default 1e-6); the old code gave a method-2 element full gravity on
  the reduced stiffness. Between start and end: factor ramps, stiffness 1
  (manual text; a stiffness ramp matched `vely` above the front but made
  `veliy` 9 % worse).
- `materi.cc` `materi()`: the stiffness factor now scales the material
  tangent `ddsdde` and the stress increment right after `set_stress()`
  (a material with the reduced stiffness). The old post-hoc scaling of
  the whole `element_matrix` / `element_lhside` left the internal force
  and the stress-projection rows unscaled: Newton updates 1e6 too large
  (Bi-CG divergence) or full-stiffness stresses stored in the soft
  elements (sigyy 349 at the top of `mesh_activate1`).

## Hardcoded parameters (F4 T8)

- Stiffness factor default 1e-6 (manual). Activation by element extreme
  vertical coordinates (not per integration point).

## Pending refactorings (F4 T8)

- `-method1` element removal; per-integration-point activation (would
  close the 2e-9 residual); `materi.cc` `strain_settlement` code still
  reads `MESH_ACTIVATE_GRAVITY_TIME` index 0 only.

