# control_mesh_remove

## Description

Deletes mesh elements according to a method. Used to remove parts of a mesh
that are completely covered by (inside) elements of other groups, e.g. to
remove embedded/surrounded volumes before an analysis.

## Usage

```
control_mesh_remove <index> -method1 element_group_0 element_group_1 [element_group_2 ...]
control_mesh_remove <index> -method3
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `method` | `-method1`: remove the elements of `element_group_0` that are completely inside elements of the other listed groups. `-method3`: remove the elements whose nodes all have an mpc. |
| `element_group_0` | `-method1`: group whose elements are candidates for removal. |
| `element_group_n` | `-method1`: groups whose elements are the reference (the containers). |

### -method1: completely located inside (manual 6.229)

An element of `element_group_0` is removed when every node of it (the mid
nodes of a quadratic element included) lies inside an active element of
the groups `element_group_1`, `element_group_2`, ..., by position: the
meshes need not share nodes (F5 T4 U9f; before, the GNU compared node
numbers, BACKLOG C126). As measured on the Professional 02-08-2026:

- the integration points of the element are tested as well (owner
  decision 2026-11-13; the Professional tests the nodes only): an element
  is removed only when its nodes and its integration points all lie inside
  reference elements, with the same tolerance. So an element over a gap
  between two reference elements, or over the outside of a reentrant
  corner, is kept when an integration point falls there
  (`control_mesh_remove_method1_gap2`, `_lcorner1`);
- this is a sampling: a gap or corner between the integration points is
  not seen (`control_mesh_remove_method1_gap1`: a gap 1 < x < 1.5 under an
  element with integration points at x = 0.42 and 1.58 is removed). A
  curved quadratic element that bulges out of a convex reference element
  between its nodes is not detected either: its integration points stay
  inside (no valid quad9 layout with all nodes in a square was found with
  an integration point outside);
- each point may lie in another reference element (an element over two
  adjacent reference elements is removed, `_union1`);
- elements without an integration rule (trusses, beams, springs) and
  interface elements are tested on their nodes only;
- a node on the boundary of a reference element is inside; the tolerance
  is about 1e-4 of the size of the reference element: 2.02e-4 in the
  [-1, 1] coordinates of a bar / quadrilateral / hexahedron (Professional:
  inside at 2e-4, outside at 2.05e-4), 0.85e-4 in the barycentric
  coordinates of a triangle / tetrahedron / prism (Professional tria3:
  inside at 0.825e-4, outside at 0.9e-4; tetrahedron and prism not
  measured);
- elements of `element_group_0` are no reference elements, also when that
  group is listed again (`-method1 2 2` removes nothing);
- 2D and 3D, linear and quadratic elements; the points are located with
  the `post_point` locator (`point_el`).

### -method3: all nodes with an mpc (manual 6.229)

An element is removed when each of its nodes has an mpc, that is, is the
constrained node (`node_0`) of an `mpc_node_number` record (the input ones
and those the GNU generates, for example from `mpc_element_group` or `mpc_linear_quadratic`,
...). A node that is only a master node has no mpc (Professional
02-08-2026). `-method3` takes no element groups: with groups the run
stops with `Error detected for data item : control_mesh_remove`.

Both methods combine with `control_mesh_remove_geometry`,
`control_mesh_remove_really -no` and `_keep_geometry` (an inactive
element stays removed while all its nodes have an mpc).

Tests in `validation-suite/test-2026` (= Professional 02-08-2026 unless a
difference below says otherwise): `control_mesh_remove_method1_inside1`,
`_boundary1`, `_outside1`, `_union1`, `_gap1`, `_gap2`, `_lcorner1`, `_group0`, `_3d1`, `_3d2`,
`_tet1`, `_quad9_1`, `_quad9_2`, `_curved1`, `_react1`, `_react2`;
`control_mesh_remove_method3_1`, `_2`, `_error1`, `_react1`, `_keep1`;
node-sharing meshes `remove_method1_1`, `remove_geometry1`, `_2`.

## Restricting the remove to a geometry

```
control_mesh_remove_geometry <index> geometry_item_name geometry_item_index
```

`control_mesh_remove_geometry` (Professional 6.230, F5 T4 U5c) restricts
the `control_mesh_remove` of the same index: an element is removed only
when every node of it lies in the geometry (as the Professional
02-08-2026). Tests `validation-suite/test-2026/remove_geometry1.dat`
(the element is inside: removed), `remove_geometry2.dat` (only two nodes
on the line: kept) and `remove_method1_1.dat` (no geometry).

## Inactive instead of deleted elements

`control_mesh_remove_really index -no` (6.232) makes the removed elements
inactive instead of deleting them; they become active again when a later
`control_mesh_remove` no longer removes them. `control_mesh_remove_keep_
geometry` (6.231), `control_mesh_remove_really_activate_factor` (6.234)
and `_activate_all` (6.233) complete it: see
[control_mesh_remove_really](control_mesh_remove_really.md) (F5 T4 U9e).
Without `control_mesh_remove_really` the removed elements are deleted, as
before (`control_mesh_remove_default1`).

The remove is evaluated on every step of a control index that has the
record (also an index without `control_timestep`).

## Example

Remove the elements of group 0 that are inside elements of group 1:

```
control_mesh_remove 0 -method1 0 1
```

After this record, every element of group 0 whose nodes all lie inside
elements of group 1 is removed from the mesh (with all its associated
element data).

## Differences with the Professional version

- A curved quadratic reference element (a mid node off the straight edge):
  the GNU tests the curved element; the Professional tests the straight
  element of its corner nodes, so a node between the straight and the
  curved edge is inside in the GNU only (`control_mesh_remove_method1_
  curved1`; DIVERGENCES REMOVE-METHOD1-CURVED).
- The GNU also tests the integration points; the Professional tests the
  nodes only and removes the elements of `control_mesh_remove_method1_gap2`
  and `_lcorner1` (DIVERGENCES REMOVE-METHOD1-INTPNT).
- The tolerances above are values inside the measured Professional
  windows, not its exact limits.
