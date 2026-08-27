# control_print_dof_line_coordinates

## Description

`control_print_dof_line_coordinates` gives the vertices of the polyline
used by [`control_print_dof_line`](control_print_dof_line.md) with the
same index (manual Professional 6.274).

The start point of the first line segment is `x_0 y_0 z_0`, the end point
of the first segment (and start of the second) is `x_1 y_1 z_1`, the end
point of the second segment is `x_2 y_2 z_2`, etc. In 1D only the
x-coordinates are needed.

## Uso

```
control_print_dof_line_coordinates 30  0.5 0. 0.5 1.
```

A polyline with two segments (a V):

```
control_print_dof_line_coordinates 33  0. 0. 1. 1. 2. 0.
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `30`      | Index of the control record. Must match the `control_print_dof_line` record. |
| `x_i y_i z_i` | Coordinates of vertex `i`. Variable length: `2*ndim` values for a single segment, `nvertices*ndim` for a polyline. At least 2 vertices are required. |

## Output

No output by itself; it defines the geometry of the line printed by
[`control_print_dof_line`](control_print_dof_line.md).

## Example

```
control_print_dof_line            30  -separate_index
control_print_dof_line_coordinates 30  0.5 0. 0.5 1.
control_print_dof_line_n          30  3
```
