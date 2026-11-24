# control_mesh_delete_element

## Description

Deletes the elements with the given numbers from the mesh. Useful to create
holes in a model or to remove parts that are not needed for the analysis.

## Usage

```
control_mesh_delete_element <index> elem_0 [elem_1 ...]
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `elem_n` | Number of the element to delete. |

## Example

Delete elements 2 and 5:

```
control_mesh_delete_element 0 2 5
```

After this record the elements 2 and 5 (with all their associated element data)
are removed from the mesh.
