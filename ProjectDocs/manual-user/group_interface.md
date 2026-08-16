# group_interface

## Description

`group_interface` marks an element group as an **interface element**
group. Interface elements model joints or discontinuities between blocks
of material (e.g. between a pile and the soil). Their strains are the
displacement differences between the two opposite sides of the element,
not field gradients.

This is the first phase of the interface family (Carril A). Currently the
elastic interface law and the Fase 3 constitutive features are
implemented:

- `group_interface_materi_elasti_stiffness kn kt,first kt,second`:
  `stress_normal = kn * strain_normal`,
  `stress_shear = kt * 2 * strain_shear`.
- `group_interface_gap gap`: physical gap of the interface. Sign
  convention: **compression = positive normal strain**. The interface is
  CLOSED (full stiffness) when the accumulated normal strain
  `strain_normal > gap`, OPEN (residual stiffness) when
  `strain_normal <= gap`. A **negative** gap is a real gap: the interface
  stays open until compression exceeds |gap|. Without the record the
  interface is always closed.
- `group_interface_materi_residual_stiffness factor`: stiffness fraction
  used when the interface is open (default 0.01).
- `group_interface_materi_plasti_tension_direct tension_limit`: tensile
  limit on the TOTAL normal force `|kn*strain_normal|`; the interface opens
  in traction when the limit is exceeded (and it was still closed).
- `group_interface_materi_plasti_mohr_coul_direct phi c phi_flow`:
  cumulative Mohr-Coulomb friction. The **presence of the record activates
  the law**: max friction = `max(c + kn*strain_normal*tan(phi), 0)`, applied
  to the TOTAL tangential force (the interface slides when the limit is
  exceeded and its tangential stiffness drops to 0 while sliding). With
  `phi=0, c=0` the limit is 0 → free sliding; without the record the
  interface is purely elastic. `phi_flow` is the dilatancy angle
  (non-associated flow): plastic slip opens the interface in both sliding
  directions.
- `group_interface_materi_memory memory_type`: memory model of the
  constitutive law. `-updated_linear` (default) recomputes the interface
  normal/tangent from the current (deformed) configuration each step;
  `-total_linear` uses the time-0 reference geometry
  (`NODE_START_REFINED`) and keeps the original interface orientation.

In 2D the interface element is a quadrilateral with 4 nodes: nodes 0-1
form side 1 and nodes 2-3 form side 2.

## Uso

Place it in the data part, in the element group definition:

```
group_type 10  -materi
group_interface 10  -yes
group_interface_materi_elasti_stiffness 10  1000.0  0.0  0.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `10`      | Element group index. |
| `-yes`    | Activate interface behaviour for this group. |

## Parámetros (Fase 3)

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_interface_gap` | `gap` | Physical gap (negative = real gap, closes under compression). Closed when strain > gap, open (residual) when strain <= gap. Default without record: always closed. |
| `group_interface_materi_residual_stiffness` | `factor` | Stiffness fraction of an open interface (default 0.01). |
| `group_interface_materi_plasti_tension_direct` | `tension_limit` | Opens in traction when the total normal force `\|kn*strain_normal\|` exceeds the limit (and the interface was closed). |
| `group_interface_materi_plasti_mohr_coul_direct` | `phi c phi_flow` | phi = friction angle (rad), c = cohesion, phi_flow = dilatancy angle (rad). The record's presence activates the cumulative Mohr-Coulomb law; phi=0,c=0 gives free sliding. |
| `group_interface_materi_memory` | `memory_type` | `-updated_linear` (default) or `-total_linear`. Memory model of the interface law; `-total_linear` fixes the normal/tangent to the time-0 geometry. |

## Related

- `group_interface_materi_elasti_stiffness index kn kt,first kt,second` —
  elastic interface stiffness (normal kn, tangential kt).

## Estado de implementación

- **Implementado (Fase 1)**: elastic interface element (2D quadrilateral),
  `group_interface`, `group_interface_materi_elasti_stiffness`.
- **Implementado (Fase 3)**: `group_interface_gap`,
  `group_interface_materi_residual_stiffness`,
  `group_interface_materi_plasti_tension_direct`,
  `group_interface_materi_plasti_mohr_coul_direct` (Mohr-Coulomb acumulativo
  con history `element_interface_force_tang`; `phi_flow` = dilatancia),
  `group_interface_materi_memory` (`-updated_linear`/`-total_linear`).
  Validados con la familia `iface_mc` (13º test de `build_safe.sh`,
  7 runs): fricción alta sostiene la carga tangencial, fricción nula
  desliza libre, tracción abre la interfaz, gap cierra bajo compresión.
- **Implementado (Fase 2)**: `control_mesh_convert` — conversion
  automatica de `-bar2` a `-quad4` para interfaces: crea los 2 nodos del
  lado opuesto de la interfaz y reconecta los elementos vecinos del otro
  lado (con `control_mesh_convert_element_group` para los grupos a un
  lado).
- **Pendiente**: `group_interface_materi_memory`, conductivity/groundflow,
  y el post-proceso de interfaz. See `ProjectDocs/DESIGN-INTERFACES.md`.
