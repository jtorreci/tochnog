# control_mesh_keep_element_group

## Description

Deletes all elements that do NOT belong to the listed groups. Useful to
visualize the model by material group.

## Usage

```
control_mesh_keep_element_group <index> group_0 [group_1 ...]
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `group_n` | Group number of the elements to keep. |

## Example

Keep only the elements of group 1, deleting all the others:

```
control_mesh_keep_element_group 0 1
```
