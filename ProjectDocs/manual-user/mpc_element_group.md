# mpc_element_group / mpc_element_group_dof / mpc_element_group_geometry (+ mpc_apply / control_mpc_apply)

## Description

Automatic multi point constraints between meshes that share a surface
but are not connected by shared nodes (manual Professional 6.860/6.864/
6.866):

```
mpc_element_group              <index> <group_0> <group_1> [<group_2> ...]
mpc_element_group_dof          <index> -velx -vely ...
mpc_element_group_geometry     <index> -geometry_line <n> ...
```

- Every node of an element of group `group_0` that is located inside an
  element of group `group_1` is TIED to that element with multi point
  constraints consistent with the shape functions at the node's
  isoparametric location in the `group_1` element (a node sitting at a
  corner of the master element gets a single master with factor 1; a
  mid-edge node gets two masters with factor 0.5, etc.).
- `-all` can be used for `group_1` to select all element groups except
  `group_0`.
- `mpc_element_group_dof` lists the dofs that are set equal (default
  without the record: every principal dof).
- `mpc_element_group_geometry` restricts the nodes of `group_0` to
  those lying on the listed geometry entities (e.g. the interface line
  of a sheet-pile wall: only the wall nodes get tied).

Generated ties use the same records as
[mpc_node_number](mpc_node_number.md) (slave dof = linear combination of
master dofs, slave treated as a known/bounded quantity) and are
re-created automatically when the mesh changes.

## Switches

- `mpc_apply <switch>` (manual 6.859, no index): global switch; with
  `-no` no mpc conditions are applied. Default `-yes`.
- `control_mpc_apply <index> <switch>` (manual 6.255): per-timestep
  override for the timestep records with the same control index.

## Example (cut-off wall under a dam, `ground19_water_under_dam`)

Two meshes meet along a vertical wall at x = 50 (element groups 2 and
3); the wall geometry keeps the meshes unmerged along it:

```
control_mesh_merge              60 -yes
control_mesh_merge_geometry_not 60 -wall
mpc_element_group               0  2 3
mpc_element_group_dof           0  -velx -vely
mpc_element_group_geometry      0  -wall
```

The wall nodes of group 2 follow the group-3 nodes exactly, so the
cut-off wall is impermeable for the mechanical field (the pore-pressure
field is already decoupled by the unmerged meshes).

## Notes

- Nodes that are ALREADY members of the node list of the containing
  master element (i.e. where the meshes share nodes) are NOT tied: the
  constraint would reduce to the node itself and freeze its equation
  (the `mpc_element_group_always -no` behaviour, manual 6.861).
- Verification: `mpc7` (two small quad4 tied into one large quad4, unit
  stress patch test) reaches `sigxx = 1.0` within `1.e-8`.
