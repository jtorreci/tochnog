# control_mesh_cut_geometry

## Description

`control_mesh_cut_geometry` cuts away the part of the mesh that lies
inside a geometry and substitutes the removed part by its nodal forces:
the elements whose nodes all lie within the geometry are deleted (with
the nodes that end up attached to no remaining element) and the force
the removed elements exerted on the surviving boundary nodes is stored
as a permanent external nodal load (`node_force`). The remaining mesh
keeps the exact equilibrium it had at the moment of the cut — the
standard excavation technique (e.g. an open-pit or a tunnel: first
compute the state of the full model, then remove the excavated part and
let the equivalent nodal forces of the removed material act on the new
boundary).

The companion record `control_mesh_cut_node_force` (the Professional
manual writes `control_mesh_cut_force`; both names are accepted) selects
per space direction whether the nodal force of the cut-away mesh is
applied to the remaining mesh (`-yes`, the default) or not (`-no`).
In 2D write two switches, in 3D three.

See `mesh_cut_1.dat` and `mesh_cut_2.dat` in the test distribution, and
`earthquake_2.dat` (dynamic analyses: cutting away the parts of the
model that are no longer needed saves computing time in calculations
with many timesteps).

## Usage

```
control_mesh_cut_geometry <index> -<geometry_item> <geometry_index>
control_mesh_cut_node_force <index> [-yes | -no] [-yes | -no] [-yes | -no]
```

The geometry item is the entity whose interior is removed (points,
lines, circles, quadrilaterals, bricks, ... — the same geometry records
used by `control_mesh_delete_geometry`). Named geometry items defined
with `start_define`/`end_define` expand to the same syntax.

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index. The cut is performed when this index becomes active (between the timestep blocks, or on the timestep with the same index). |
| `-<geometry_item> <geometry_index>` | The geometry that bounds the cut-away part. |
| switches | `control_mesh_cut_node_force`: per direction, whether the nodal force of the cut-away mesh is applied to the remaining mesh (`-yes`) or not (`-no`). Missing record: all directions applied. |

## Example

A bar fixed at the left edge is elongated by pulling the right edge;
then the right half is cut away and the left half must remain in
equilibrium (1D, `mesh_cut_1.dat`):

```
bounda_dof    10  -left_edge -disx
bounda_dof    20  -right_edge -disx

control_timestep              20  1.e-1 2

control_mesh_cut_geometry     30  -right_half
control_mesh_cut_node_force   30  -yes

control_timestep              40  1.e-1 1
```

After the cut (control index 30, at the end of the loading block) the
elements of the right half are removed and the substituted nodal force
keeps the remaining bar stretched: the stress at a post point of the
left half stays at its pre-cut value.

## Differences with the Professional version

- The nodal-force substitution is computed from the element
  integration-point stresses (`ELEMENT_DOF`), which requires
  `materi_stress` in the `initia` section and the default
  `options_element_dof -yes`. The supported elements are the
  tensor-product isoparametric family: `-bar2`/`-bar3` (1D),
  `-quad4`/`-quad9` (2D), `-hex8`/`-hex27` (3D). Other element types in
  the cut geometry abort with a clear message (no deletion without the
  equilibrium substitution). Axisymmetric groups are not supported by
  the cut integration.
- The Professional stores the substituted `node_force` with the sign of
  the external load it applies; this GNU version stores the negated
  value (its `node_force` consumption applies a stored `+F` as a `-x`
  load — see the developer manual; the equilibrium after the cut is
  identical).
