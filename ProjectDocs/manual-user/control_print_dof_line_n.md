# control_print_dof_line_n

## Description

`control_print_dof_line_n` determines how many points are printed along
the line of [`control_print_dof_line`](control_print_dof_line.md)
(manual Professional 6.279).

The points are distributed with equal spacing over the TOTAL length of
the polyline (including the start and end points). For `n = 5` on a
single segment you get the points at fractions 0, 1/4, 1/2, 3/4 and 1 of
the segment.

Default, when the record is not given: `5` (the legacy GNU
`post_line_n` default; the Professional manual does not state a
default).

## Uso

```
control_print_dof_line_n 30  5
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `30`      | Index of the control record. Must match the `control_print_dof_line` record. |
| `n`       | Number of points printed along the line, `n >= 1`. `n = 1` prints only the start point. |

## Output

No output by itself; it defines the number of data lines written per
call by [`control_print_dof_line`](control_print_dof_line.md).

## Example

```
control_print_dof_line           30  -separate_index
control_print_dof_line_coordinates 30  0.5 0. 0.5 1.
control_print_dof_line_n         30  3
```
