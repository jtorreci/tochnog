# control_mesh_change_element_group

## Description

Changes the group of all elements belonging to `group_from` to `group_to`.
Requires that the elements have `element_group` defined explicitly.

## Usage

```
control_mesh_change_element_group <index> group_from group_to
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index. Must NOT be equal to a timestep index. |
| `group_from` | Group number of the elements to change. |
| `group_to` | New group number assigned to those elements. |

## Example

Move all elements of group 0 to group 1:

```
control_mesh_change_element_group 0 0 1
```
