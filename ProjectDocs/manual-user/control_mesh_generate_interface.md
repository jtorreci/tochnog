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

## Elementos generados

| Shared face | Generated element |
|-------------|-------------------|
| 2D, 2 nodes (between two `quad4`) | `-quad4` |
| 3D, 3 nodes (between two `tet4`) | `-prism6` |
| 3D, 4 nodes (between two `hex8`) | `-hex8` |

## Related

- `group_interface` — the element group of the generated interface
  elements (must be set on `eg_i`).
- `control_mesh_convert` — converts low-dimensional interface elements to
  their isoparametric equivalent.
- `control_mesh_generate_interface_geometry index geometry_item_name
  geometry_item_index` — restrict generation to a geometry.

## Estado de implementación

- **Implementado**: `control_mesh_generate_interface` (2D `quad4`, 3D
  `prism6`/`hex8`) + `control_mesh_generate_interface_geometry`.
  Validado con `iface_gen` (dos quad4 con nodos duplicados en x=1 →
  interfaz `{2 4 5 7}` generada; la fricción MC sostiene la carga
  tangencial), `iface_gen_geom` (geometría que cubre la interfaz → se
  genera) e `iface_gen_geom_off` (geometría desplazada → no se genera,
  el bloque se mueve libre).
- **Pendiente**: `control_mesh_generate_interface_method`
  (`-element_geometry`); elementos cuadráticos (quad6, tria12, quad18).
