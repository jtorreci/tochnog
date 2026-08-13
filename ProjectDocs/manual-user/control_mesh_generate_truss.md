# control_mesh_generate_truss

## Description

`control_mesh_generate_truss` generates truss elements between nodes that
are neighbours in space (that is, nodes connected by an isoparametric
finite element). Only nodes located on the specified geometry entity are
used. This is convenient for adding trusses to a mesh that already has
isoparametric elements.

The generated trusses get an `element_group` record with the given value,
so you can attach truss properties (`group_truss_elasti_young`, etc.) to
that group.

## Uso

Place it in the data part:

```
control_mesh_generate_truss 10  1 -geometry_line 1
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `10`      | Index of the control record. |
| `1`       | Element group assigned to the generated trusses. |
| `-geometry_line 1` | Geometry entity; only nodes on it are used. |

## Related

- `control_mesh_generate_beam` — same, but generates beam elements.
- `control_mesh_generate_trussbeam` — same, but generates truss-beam
  elements.
- `control_mesh_generate_truss_beam_loose` — with `-yes`, the generated
  elements are not connected to the existing nodes (new nodes are
  generated).
- `control_mesh_generate_truss_beam_macro` — truss/beam macro.

## Example

A 1D bar with a truss generated between its two nodes (the right node gets
a force leading to a displacement of 0.5):

```
geometry_line 1 0. 1. 1.e-4
node 1  0
node 2  1
element 1  -bar2 1 2
control_mesh_generate_truss 10  1 -geometry_line 1
control_timestep            20  1. 1.
target_item 0   -node_dof 2 -disx
target_value 0   0.5 1.e-2
```
