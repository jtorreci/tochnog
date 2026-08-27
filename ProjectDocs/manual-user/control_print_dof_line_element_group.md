# control_print_dof_line_element_group

## Description

`control_print_dof_line_element_group` restricts the interpolation of
[`control_print_dof_line`](control_print_dof_line.md) to the given
element groups (manual Professional 6.275): each point of the line is
only searched in elements of those groups.

A point that lies on an element OUTSIDE the listed groups is not printed.

## Uso

```
control_print_dof_line_element_group 34  1
```

Only elements of group 1 are searched; a point of the line that only
lies in an element of another group is omitted.

Several groups can be listed:

```
control_print_dof_line_element_group 34  0 1 3
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `34`      | Index of the control record. Must match the `control_print_dof_line` record. |
| `element_group_i` | Element groups (integers) in which the line points are searched. Variable length. |

## Output

No output by itself; it restricts the interpolation of
[`control_print_dof_line`](control_print_dof_line.md).

## Example

```
control_print_dof_line               34  -separate_index
control_print_dof_line_coordinates   34  0. 0.5 2. 0.5
control_print_dof_line_n             34  3
control_print_dof_line_element_group 34  1
```
