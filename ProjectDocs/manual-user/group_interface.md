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
- `group_interface_gap gap`: initial empty space between the sides; the
  interface only generates stresses when closed (normal strain < gap).
- `group_interface_materi_residual_stiffness factor`: stiffness fraction
  used when the interface is open (default 0.01).
- `group_interface_materi_plasti_tension_direct tension_limit`: tensile
  limit; beyond it the interface opens.
- `group_interface_materi_plasti_mohr_coul_direct phi c phi_flow`:
  Mohr-Coulomb friction (max friction = c + Fn*tan(phi)).

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

## Related

- `group_interface_materi_elasti_stiffness index kn kt,first kt,second` —
  elastic interface stiffness (normal kn, tangential kt).

## Estado de implementación

- **Implementado (Fase 1)**: elastic interface element (2D quadrilateral),
  `group_interface`, `group_interface_materi_elasti_stiffness`.
- **Implementado (Fase 3)**: `group_interface_gap`,
  `group_interface_materi_residual_stiffness`,
  `group_interface_materi_plasti_tension_direct`,
  `group_interface_materi_plasti_mohr_coul_direct`.
  Validación completa de gap; Mohr-Coulomb/tension en validación de humo.
- **Pendiente**: `control_mesh_convert` (automatic conversion of bar2 ->
  quad4 etc.), `group_interface_materi_memory`, conductivity/groundflow,
  and the interface post-processing. See `ProjectDocs/DESIGN-INTERFACES.md`.
