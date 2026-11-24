# control_mesh_keep_element

## Description

Deletes all elements of the mesh EXCEPT the ones listed. Useful to visualize
only a specific set of elements in isolation.

## Usage

```
control_mesh_keep_element <index> elem_0 [elem_1 ...]
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `elem_n` | Number of the element to keep. |

## Example

Keep only element 2, deleting all the others:

```
control_mesh_keep_element 0 2
```
