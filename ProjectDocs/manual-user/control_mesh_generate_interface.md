# control_mesh_generate_interface

## Description

`control_mesh_generate_interface` generates **interface elements**
between two element groups that share a common face. It is the mesh
generation counterpart of `group_interface`: it creates the interface
elements that connect two blocks of material (e.g. a pile and the soil,
or two meshes with duplicated nodes on the interface face).

The two blocks must have **duplicated (spatially coincident) nodes** on
the interface face. For each element pair (one in each group) that shares
a complete face, an interface element is generated with the nodes of BOTH
sides and assigned to the given interface element group.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_timestep` record:

```
control_mesh_generate_interface 0  10 0 0
```

The record is a list of triples `eg_i eg_a eg_b`:
- `eg_i`: element group of the generated interface element (must have
  `group_interface -yes`).
- `eg_a`: element group on one side of the interface.
- `eg_b`: element group on the opposite side.

Up to several interfaces can be generated in one record:
`eg_0 eg_00 eg_01 eg_1 eg_10 eg_11 ...`. For connected interfaces, put
the touching pairs in ONE record.

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `0`       | Index of the control record. Must match the `control_timestep` index. |
| `10 0 0`  | Triple: interface group `10`, side-A group `0`, side-B group `0`. |

## Restricción por geometría

`control_mesh_generate_interface_geometry` restricts the generation to a
geometry:

```
geometry_quadrilateral 1  0.9 0.  1.1 0.  1.1 1.  0.9 1.  0.01
control_mesh_generate_interface_geometry 0  -geometry_quadrilateral 1
```

Only element pairs whose shared nodes lie inside the geometry are
converted.

## Selección y generación por `element_geometry`

`control_mesh_generate_interface_method index method_select method_generate`
changes how elements are selected and how the generated interface is
assigned:

- `method_select = -element_geometry`: the pairs `eg_a eg_b` select
  elements by `element_geometry` instead of `element_group`.
- `method_generate = -element_geometry`: the generated interface element
  receives an `element_geometry` record (instead of `element_group`).

`element_geometry index geometry_set` assigns a geometrical set number to
an element:

```
element_geometry 1  20
element_geometry 2  30
control_mesh_generate_interface 0  10 20 30
control_mesh_generate_interface_method 0  -element_geometry -element_geometry
```

## Elementos generados

| Shared face | Generated element |
|-------------|-------------------|
| 2D, 2 nodes (between two `quad4`) | `-quad4` |
| 3D, 3 nodes (between two `tet4`) | `-prism6` |
| 3D, 4 nodes (between two `hex8`) | `-hex8` |

**Elementos cuadráticos**: when the contacting elements are quadratic
(`quad9`, `tet10`, `hex27`, `bar3`), the shared face has mid-side nodes.
The face is **subdivided into linear interface elements** so all face
nodes (including mid-side and centre) are coupled:

| Shared face | Sub-division |
|-------------|--------------|
| 2D quadratic edge (3 nodes, `quad9`/`bar3`) | 2 `-quad4` |
| 3D tria6 (6 nodes, `tet10`) | 4 `-prism6` |
| 3D quad9 (9 nodes, `hex27`) | 4 `-hex8` |

This is the current behaviour. Native quadratic interface elements
(`quad6`/`tria12`/`quad18`) are a documented future option.

## Related

- `group_interface` — the element group of the generated interface
  elements (must be set on `eg_i`).
- `control_mesh_convert` — converts low-dimensional interface elements to
  their isoparametric equivalent.
- `control_mesh_generate_interface_geometry index geometry_item_name
  geometry_item_index` — restrict generation to a geometry.
- `control_mesh_generate_interface_method index method_select
  method_generate` — select by `-element_geometry` and/or generate an
  `element_geometry` record.
- `element_geometry index geometry_set` — assign a geometrical set to an
  element.

## Estado de implementación

- **Implementado**: `control_mesh_generate_interface` (2D `quad4`, 3D
  `prism6`/`hex8`) + `control_mesh_generate_interface_geometry` +
  `control_mesh_generate_interface_method` (`-element_geometry`) +
  **subdivisión de caras cuadráticas** (quad9/tet10/hex27/bar3).
  Validado con `iface_gen` (dos quad4 con nodos duplicados en x=1 →
  interfaz `{2 4 5 7}` generada; la fricción MC sostiene la carga
  tangencial), `iface_gen_geom` (geometría que cubre la interfaz → se
  genera), `iface_gen_geom_off` (geometría desplazada → no se genera, el
  bloque se mueve libre), `iface_gen_method` (selección por
  `element_geometry`), `iface_gen_method_gen` (generación con
  `element_geometry`) e `iface_gen_quad9` (dos `quad9` con arista
  cuadrática → 2 interfaces `quad4`; la fricción MC sostiene).
- **Pendiente / futuro**: interfaces cuadráticas nativas
  (`quad6`/`tria12`/`quad18`); soporte serendipito (`hex20`/`quad8`).
