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

## Example

Remove the elements of group 0 that are inside elements of group 1:

```
control_mesh_remove 0 -method1 0 1
```

After this record, every element of group 0 whose nodes are all contained in
an element of group 1 is removed from the mesh (with all its associated element
data).
