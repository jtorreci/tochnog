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
- Follows the `spring.cc` pattern: assembly on `vel_indx` (velocity dof).
  The nodal force is `-sign*stress*dir` (principle of virtual work), the
  matrix/lhside from `K*dtime`.
- **Sign convention (2026-08-13)**: the nodal force must be
  `-sign*(stress*dir)`; with the opposite sign the interface pushes
  instead of resisting (kn high increased the displacement).

## Validación

- Test 2 bloques (`/tmp/iface2.dat`): un bloque izquierdo fijo y un
  bloque derecho empujado, conectados por la interfaz. Verificado:
  - kn=100 → velix(nodo6)=0.044
  - kn=1000 → 0.0018
  - kn=1e6 → -0.0032
  - kn=0.001 → ≈1.0 (deslizamiento libre)
  - Sin interfaz (bloques soldados) → -0.0066
  El límite kn→∞ tiende al modelo soldado, kn→0 al deslizamiento libre —
  comportamiento físico correcto.
- La fuerza usa la velocidad relativa entre lados `(v_side2-v_side1)*dtime`
  (incremento de desplazamiento), no `dis_indx` (que es -1 con
  velocity_integrated).

## Pendiente / validación

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
