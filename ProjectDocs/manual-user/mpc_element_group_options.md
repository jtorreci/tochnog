# mpc_element_group options

`mpc_element_group_always`, `mpc_element_group_closest`,
`mpc_element_group_coord_geometry`, `mpc_element_group_eps_iso`,
`mpc_element_group_keep` and `control_mpc_element_group`.

## Description

`mpc_element_group index element_group_0 element_group_1` (manual
Professional 6.858, see [mpc_element_group](mpc_element_group.md)) ties every
node of an element of `element_group_0` that lies inside an element of
`element_group_1` to that element with multi point constraints. These
records control which nodes get a tie and when the ties are created
(manual Professional 6.256, 6.859-6.865). They are useful to connect
meshes that were made apart, for example soil nodes to a pile: tie only
the soil nodes beside the pile shaft, and only above its toe.

| Record | Meaning | Default |
|--------|---------|---------|
| `mpc_element_group_always index switch` (6.859) | `-no`: a node that is a member of the node list of an element of `element_group_1` (the two groups are connected there) gets no tie. `-yes`: such a node is still tied to another `element_group_1` element that contains it. A master element never ties a node to itself. | `-yes` |
| `mpc_element_group_closest index switch` (6.860) | `-yes`: a node that lies in no element of `element_group_1` is tied to the closest one, with the shape functions at the closest location in that element. | `-no` |
| `mpc_element_group_coord_geometry index switch` (6.861) | `-yes`: only nodes whose last coordinate (y in 2D, z in 3D) is not below the lowest node of `element_group_1` are tied. | `-no` |
| `mpc_element_group_eps_iso index eps` (6.863) | Tolerance on the isoparametric coordinates (range -1 to 1) of the `element_group_1` element below which a node counts as inside it. | `1.e-4` |
| `mpc_element_group_keep index switch` (6.865) | `-yes`: the ties created for this index stay for the rest of the calculation; they are not created again after a mesh change or by `control_mpc_element_group`. Ties whose nodes are deleted are dropped. | `-no` |
| `control_mpc_element_group index switch` (6.256) | `index` is a control index. `-yes`: the ties of all `mpc_element_group` records are created again at every timestep of that control index (for example after a `change_dataitem` of the geometry of `mpc_element_group_geometry`). `-no`: only when the mesh has changed. | `-no` |

The `index` of the first five records is the index of the
`mpc_element_group` record they belong to; that record must exist.

## Usage

```
mpc_element_group                10 soil_group pile_group
mpc_element_group_geometry       10 -geometry_line 1
mpc_element_group_closest        10 -yes
mpc_element_group_coord_geometry 10 -yes
mpc_element_group_dof            10 -velx -vely
```

## Example

```
node 1 0.0 0.25
node 2 0.0 1.0
node 3 1.0 0.25
node 4 1.0 1.0
element 1 -quad4 1 2 3 4
node 5 1.0 0.0
node 6 1.0 0.5
node 7 1.0 1.0
...
element_group 1 1
element_group 2 2
element_group 3 2
mpc_element_group 10 2 1
mpc_element_group_dof 10 -velx -vely
geometry_line 1 1.0 -1.0 1.0 2.0 1.e-6
mpc_element_group_geometry 10 -geometry_line 1
mpc_element_group_closest 10 -yes
mpc_element_group_coord_geometry 10 -yes
```

Node 5 (y = 0) lies below the lowest node of group 1 (y = 0.25) and gets
no tie; node 6 gets the factors 2/3 1/3 on nodes 3 and 4, node 7 the
factor 1 on node 4. Tests: `validation-suite/test-2026/mpc_element_group_*.dat`
and `control_mpc_element_group_*.dat`.

## Differences with the Professional version

Compared with the Professional 02-08-2026 binary:

- `mpc_element_group_always -no`: the Professional ignores it and also ties
  nodes that are members of an `element_group_1` element. The GNU follows
  the manual.
- `mpc_element_group_coord_geometry -yes`: the Professional takes the lowest
  node only among the `element_group_1` nodes on
  `mpc_element_group_geometry`, and ignores the switch without that record.
  The GNU follows the manual (lowest node of `element_group_1`).
- `mpc_element_group_eps_iso`: without the record the Professional accepts
  nodes up to about 1.e-3 outside a unit element (the GNU follows the manual
  default 1.e-4), and with a tolerance of 1.e-2 or more it ties no node
  outside the element at all. The GNU factors of an accepted outside node
  are extrapolated, the Professional ones clamped (BACKLOG C127).
- `mpc_element_group_keep -yes`: ignored by the Professional.
- `control_mpc_element_group -no`: the Professional creates the ties again
  at every timestep anyway. The GNU follows the manual.
- `mpc_element_group_closest`: equal on rectangular quad4 / hex8 and on
  tria3 masters; on a distorted quad4 the Professional uses an
  approximation of the closest point, and on tet4 masters it ignores the
  fourth corner.

All are listed in `ProjectDocs/DIVERGENCES-PROFESSIONAL.md`.
