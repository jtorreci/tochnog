# control_mesh_remove

## Description

Deletes mesh elements according to a method. Used to remove parts of a mesh
that are completely covered by (inside) elements of other groups, e.g. to
remove embedded/surrounded volumes before an analysis.

Partially implemented: only `-method1` works in this version.

## Usage

```
control_mesh_remove <index> <method> element_group_0 [element_group_1 ...]
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `method` | `-method1`: remove the elements of `element_group_0` that are completely inside elements of the other listed groups. |
| `element_group_0` | Group whose elements are candidates for removal. |
| `element_group_n` | Groups whose elements are used as the containment reference. |

Note: `-method3` (remove elements where all nodes have an mpc) is NOT
implemented in this version; only `-method1` is accepted.

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

After this record, every element of group 0 whose nodes are all contained in
an element of group 1 is removed from the mesh (with all its associated element
data).
