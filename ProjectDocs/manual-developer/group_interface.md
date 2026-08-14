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

## Ley constitutiva (Fase 3)

- **Sign convention (RF-5)**: compression = `strain_normal` POSITIVE,
  traction NEGATIVE. This is the convention used by gap, tension_direct
  and `max_fric` (`kn*strain_normal` grows with compression). The old doc
  said "compression negative" — wrong, corrected 2026-08-14.
- **Gap** (`group_interface_gap`): the interface is OPEN (residual
  stiffness only) when the accumulated normal strain `strain_normal <= gap`,
  CLOSED (full stiffness) when `strain_normal > gap`. Compression
  (strain > 0) always closes. A physical gap is a NEGATIVE value: the
  interface stays open until compression exceeds |gap|. The accumulated
  strain is stored in `ELEMENT_INTERFACE_STRAIN_NORMAL`. If no record is
  given the interface is always closed (default gap = -1e20; the old +1e20
  default with the inverted condition left the interface always open).
- **Residual stiffness** (`group_interface_materi_residual_stiffness`):
  fraction of the original stiffness used when the interface is open
  (default 0.01).
- **Tension limit** (`group_interface_materi_plasti_tension_direct`): the
  interface opens in TRACTION when the accumulated TOTAL normal force
  `|Fn_total| = |kn*strain_normal| > tension_limit` (requires
  `strain_normal < 0`) and only if it was still closed
  (`stiff_normal == kn`). On opening: residual stiffness and `force_norm`
  capped at `±tension_limit` (signed by the opening direction).
- **Cumulative Mohr-Coulomb** (`group_interface_materi_plasti_mohr_coul_direct
  phi c phi_flow`): active by the PRESENCE of the record (D2). With
  phi=0,c=0 the limit is 0 → free sliding; without the record the interface
  stays purely elastic (Fase 1). The limit applies to the TOTAL tangential
  force, stored in the history `ELEMENT_INTERFACE_FORCE_TANG` (spring.cc
  pattern, GET VERSION_NORMAL with GET_IF_EXISTS, PUT VERSION_NEW):
  - `trial = f_t_old + kt1*2*du_tang` (accumulated across steps; the
    step loop copies NORMAL→NEW at step start, NEW→NORMAL at step close)
  - `max_fric = max(c + kn*strain_normal*tan(phi), 0)` (floor at 0, D5)
  - clamp: `f_t = clamp(trial, ±max_fric)`; if plastified →
    `stiff_tang = 0` (consistent tangent, D3, avoids Newton oscillation
    at the elastic/plastic boundary)
  - `stress_shear = f_t - f_t_old` (the rhs carries the INCREMENT; without
    MC this equals `kt1*2*du_tang` exactly — Fase 1 backward compatible)
- **phi_flow = dilatancy (RF-4, non-associated flow)**: if plastified and
  `phi_flow > 0`, slip OPENS the interface in both sliding directions:
  `strain_normal += -|du_tang| * tan(phi_flow)` (magnitude, not signed
  du_tang — a signed flow would close one direction, anti-physical). Feeds
  back into the normal history → gap/tension/max_fric of the next step.

## Validación

- **Fase 1 elástica**: test 2 bloques (`/tmp/iface2.dat`) — kn=100→0.044,
  kn=1000→0.0018, kn=1e6→-0.003 (aprox soldado), kn=0.001→≈1.0 libre,
  sin interfaz→-0.0066. Límites físicos correctos.
- **Fase 3 validada (familia `iface_mc`, 13º test de build_safe.sh, 6 runs)**:
  - (a/a') fricción alta sostiene: phi=45°, c=0, Fy=5 < límite ≈ 80 →
    vely(nodo 6) ≈ 0 (0.20 en 10 pasos / 0.006 en 1 paso — invarianza de
    nº de pasos OK, mismos targets).
  - (b/b') fricción nula desliza: phi=0, c=0 con record PRESENTE → límite 0
    → `element_interface_force_tang`(2) = 0 (deslizamiento libre; fallback
    cinemático: vely=5 prescrita en los nodos 5-8).
  - (c) tracción abre: tension_limit=1.0, Fx=+10 → velix(nodo 6) = +3.88
    (strain_normal < 0; el código viejo exigía compresión y nunca abría).
  - (d) gap cierra con compresión: gap=0.001, Fx=-10 → strain 0.138 > gap →
    velix(nodo 6) ≈ 0 (−0.06; el código viejo abría bajo compresión).
  - Discriminadores vs el código viejo: (b)/(b') elástico puro →
    force_tang ~10⁴ ≫ 1.0 FALLA; (c) nunca abre → velix −13.66 FALLA;
    (d) abre con compresión → velix grande FALLA.

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
- `control_mesh_convert` (bar2 -> quad4) implementado en la Fase 2 (commit
  `490545b`); el quad4 de interfaz también puede definirse manualmente.
- Solo 2D; el caso 3D (normal en z, 2 tangentes) está esbozado pero no
  validado.

## Detalles

- `normal`/`tangent` se calculan de la geometría: tangent a lo largo del
  lado 1, normal perpendicular (2D).
- La deformación usa el incremento `du_new - du_old` (acumulativo).
