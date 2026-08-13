# group_interface

## Implementación

- **Elemento**: `interface_element()` in `interface.cc` (new file, added
  to the makefile as `interface.$(OBJ)`). Invoked from `elem()` in
  `elem.cc`, in the structural-elements branch, when
  `db_active_index( GROUP_INTERFACE, element_group )` is true.
- **Keywords** (data_class GROUP_TYPE) registered in `database.cc`:
  - `group_interface` (INTEGER, length 1, required GROUP_TYPE).
  - `group_interface_materi_elasti_stiffness` (DOUBLE, length 3, required
    GROUP_INTERFACE).
- **New enums**: `GROUP_INTERFACE`, `GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS`
  in `tochnog.h` / `tochnog-mod.h` (kept in sync), placed between
  `GROUP_INTEGRATION_POINTS` and `GROUP_MATERI_DAMAGE_MAZARS`.

## Modelo físico

- Interface strains = displacement differences between the two opposite
  sides (not gradients). In 2D the element is a quadrilateral with 4
  nodes; nodes {0,1} form side 1 and nodes {2,3} side 2.
- Elastic law:
  - `stress_normal = kn * strain_normal`
  - `stress_shear = kt * 2 * strain_shear`
  where strain = displacement difference / thickness (thickness absorbed
  in kn/kt, so thickness = 1).
- The displacement dof is `veli_indx` (velocity_integrated) or `dis_indx`
  (displacement), selected via `materi_velocity_integrated`.

## Ensamblaje

- The stiffness matrix is `[K -K; -K K]` on the displacement dofs of the
  two sides, with K = kn (normal) + kt*2 (tangential) projected onto the
  interface normal/tangent.
- Follows the `spring.cc` pattern: assembly on `vel_indx` (velocity dof),
  `element_rhside` from the stress, `element_matrix`/`element_lhside`
  from K*dtime.

## Pendiente / validación

- **Validación física del acoplamiento**: la Fase 1 se implementa y el
  ensamblaje se verifica (matriz y fuerza se generan), pero el
  comportamiento del test de 2 bloques con kn alta requiere afinado
  (el desplazamiento del bloque derecho no se reduce con kn alta como
  se espera). Pendiente de análisis del acoplamiento con el solver.
- El test aislado (1 elemento de interfaz con lados prescritos) no valida
  bien porque los dofs prescritos no dejan que la interfaz frene.
- `control_mesh_convert` (bar2 -> quad4 etc.) no está implementado — el
  quad4 de interfaz debe definirse manualmente.
- Solo 2D; el caso 3D (normal en z, 2 tangentes) está esbozado pero no
  validado.

## Detalles

- `normal`/`tangent` se calculan de la geometría: tangent a lo largo del
  lado 1, normal perpendicular (2D).
- La deformación usa el incremento `du_new - du_old` (acumulativo).
