# control_mesh_truss_distribute_mpc

`control_mesh_truss_distribute_mpc` and its options `_air`, `_dof`,
`_element_group_truss`, `_element_group_isoparametric`, `_exact`,
`_exact_minimal_length`, `_exact_minimal_length_connect`,
`_geometry_truss` and `_geometry_isoparametric` (manual Professional
6.245-6.254).

## Description

With `control_mesh_truss_distribute_mpc index -yes` the nodes of the
`-truss` elements are fixed with multi point constraints (mpc's) to the
isoparametric elements the trusses run through. This models reinforcement
bars in a concrete embedment: the bars follow the displacements (and the
temperature, if present) of the concrete, without a mesh that shares the
bar nodes.

Each tied node gets one `mpc_node_number` / `mpc_node_factor` record per
dof: the masters are the nodes of the isoparametric element that contains
it, the factors their shape functions at the node position. The records
are written to the `.dbs` with an `mpc_from_truss_distribute_mpc -yes`
marker and work as user `mpc_node_number` records.

The command runs once, when its control index is processed (a control
index without `control_timestep` is fine). Every distributed truss is
replaced by new elements numbered after the highest element (the old
element is deleted, its group and other element data are copied).

With `control_mesh_truss_distribute_mpc_exact index -yes` a truss is split
where it enters a next isoparametric element, so each new truss lies in one
element. This helps when the trusses are long compared with the
isoparametric elements. The new nodes at the split points are numbered
after the highest node.

| Record | Meaning | Default |
|--------|---------|---------|
| `control_mesh_truss_distribute_mpc index switch` (6.245) | `-yes`: distribute the trusses at control index `index`. Only one record is allowed. | - |
| `control_mesh_truss_distribute_mpc_exact index switch` (6.250) | `-yes`: split the trusses at the element entries. | `-no` |
| `control_mesh_truss_distribute_mpc_exact_minimal_length index length` (6.251) | Split trusses shorter than `length` are not generated. | 1.e-11 times the truss length |
| `control_mesh_truss_distribute_mpc_exact_minimal_length_connect index switch` (6.252) | `-yes`: when a short truss between two generated trusses is dropped, the next truss starts at the node where the dropped one started (connected). `-no`: it starts at a new node (loose). Only trusses made from one old truss are connected. | `-yes` |
| `control_mesh_truss_distribute_mpc_air index switch` (6.246) | `-yes`: trusses whose centre is not inside an isoparametric element are also generated. `-no`: they are not generated. | `-yes` |
| `control_mesh_truss_distribute_mpc_dof index dof_0 dof_1 ...` (6.247) | The dofs that are tied, e.g. `-velx -vely -temp`. | all principal dofs |
| `control_mesh_truss_distribute_mpc_element_group_truss index group_0 group_1 ...` (6.248) | Only the trusses of these element groups are distributed. | all trusses |
| `control_mesh_truss_distribute_mpc_geometry_truss index name_0 index_0 ...` (6.253) | Only the trusses with all nodes in one of these geometries (a `geometry_set` is allowed). | all trusses |
| `control_mesh_truss_distribute_mpc_element_group_isoparametric index group_0 ...` (6.249) | Only the isoparametric elements of these groups are used. | all |
| `control_mesh_truss_distribute_mpc_geometry_isoparametric index name_0 index_0 ...` (6.254) | Only the isoparametric elements in these geometries are used. | all |

Only one of `_element_group_truss` and `_geometry_truss` may be given, and
only one of the two `_isoparametric` records. When a truss record and an
isoparametric record are both given they need the same number of values:
the first truss group (or geometry) is distributed over the first
isoparametric group (or geometry), the second over the second, and so on.

Which nodes are tied:

- a truss node that is a node of a used isoparametric element is not tied
  (the meshes are connected there);
- a node inside a used isoparametric element is tied to the first one (in
  element order) that contains it;
- without `_exact` a node outside every used element (in the air) is not
  tied; with `_exact` the air end of a split truss whose centre lies in an
  element is tied to that element at the point where the truss leaves it.

The isoparametric elements are the continuum elements of the full
dimension: `-bar2` / `-bar3` in 1D, `-tria3` / `-tria6` / `-quad4` /
`-quad9` in 2D, `-tet4` / `-tet10` / `-hex8` / `-hex27` / `-prism6` in 3D,
outside interface groups. The entry points are computed with straight
element sides.

## Usage

```
control_mesh_truss_distribute_mpc                     10 -yes
control_mesh_truss_distribute_mpc_exact               10 -yes
control_mesh_truss_distribute_mpc_element_group_truss 10 1
control_mesh_truss_distribute_mpc_dof                 10 -velx -vely
```

## Example

```
number_of_space_dimensions 1
materi_velocity
materi_displacement
materi_stress
end_initia
node 1 0
node 2 10
node 3 20
node 4 2
node 5 18
element 1 -bar2 1 2
element 2 -bar2 2 3
element 3 -truss 4 5
element_group 3 1
group_type 1 -truss
group_truss_area 1 0.1
group_truss_elasti_young 1 1.
control_mesh_truss_distribute_mpc       10 -yes
control_mesh_truss_distribute_mpc_exact 10 -yes
end_data
```

The truss is replaced by element 4 (nodes 4 and 6) and element 5 (nodes 6
and 5); the new node 6 lies at x = 10, where the truss enters the second
bar. Node 4 gets the factors 0.8 / 0.2 on nodes 1 / 2, node 5 0.2 / 0.8 on
nodes 2 / 3 and node 6 the factor 1 on node 2. Tests:
`validation-suite/test-2026/truss_distribute_mpc_*.dat` and the copies of
the Professional examples `truss10.dat` ... `truss15.dat`.

## Differences with the Professional version

Compared with the Professional 02-08-2026 binary (DIVERGENCES-PROFESSIONAL.md):

- `_air` default: the manual says `-yes`, the Professional behaves as `-no`
  when the record is missing. The GNU follows the manual.
- `_air -no` without `_exact`: the Professional deletes a truss as soon as
  one of its nodes is outside the isoparametric elements; the GNU follows
  the manual (only when the centre of the truss is outside).
- `_air -yes` with `_exact`: the Professional also splits the truss where
  it leaves into the air and leaves the air node untied; the GNU splits
  only at element entries (manual) and ties the air node at the exit point
  when the piece centre is inside an element.
- `_element_group_isoparametric` / `_geometry_isoparametric`: the
  Professional ignores the restriction; the GNU follows the manual.
- Numbering: the Professional leaves gaps in the new node and element
  numbers (and in the mpc indices); the GNU numbers them compactly. The
  mesh, the ties and the results are the same. A target on a new element
  number (Professional `truss13`, element 8) uses the GNU number.
- Zero-factor masters are left out of the records (the Professional lists
  them); the tie is the same.
- `_dof`: the manual heading shows no index; the record takes the control
  index as in the Professional.
- Generating `-trussbeam` elements instead of trusses (manual 6.245, "as a
  special option") is not implemented: the manual names no record for it.
