# group_interface

## Description

`group_interface` marks an element group as an **interface element**
group. Interface elements model joints or discontinuities between blocks
of material (e.g. between a pile and the soil). Their strains are the
displacement differences between the two opposite sides of the element,
not field gradients.

This is the first phase of the interface family (Carril A). Currently only
the elastic interface law is implemented:

- `group_interface_materi_elasti_stiffness kn kt,first kt,second`:
  `stress_normal = kn * strain_normal`,
  `stress_shear = kt * 2 * strain_shear`.

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
- **Pendiente**: `control_mesh_convert` (automatic conversion of bar2 ->
  quad4 etc.), the plastic interface laws (`_mohr_coul_direct`,
  `_tension_direct`), gap, memory, and the interface post-processing.
  See `ProjectDocs/DESIGN-INTERFACES.md`.
