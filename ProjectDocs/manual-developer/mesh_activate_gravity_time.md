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
