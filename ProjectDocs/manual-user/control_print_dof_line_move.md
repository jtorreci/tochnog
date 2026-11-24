# control_print_dof_line_move

## Description

`control_print_dof_line_move` moves the coordinates of the line of
[`control_print_dof_line`](control_print_dof_line.md) with the same
index along with the velocity field (manual Professional 6.278). With
this option you can follow, with a `control_print_dof_line`, the dofs of
material particles: after each print the line vertices are displaced by
the interpolated velocity times the time step.

Use this option only when `materi_velocity` is initialised in the
initialisation part (otherwise the line is not moved).

## Uso

```
control_print_dof_line_move 40  -yes
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `40`      | Index of the control record. Must match the `control_print_dof_line` record. |
| switch    | `-yes`: move the line coordinates with the velocity field. `-no` (default): keep them fixed. |

## Output

No output by itself; it updates the stored line coordinates after every
print of [`control_print_dof_line`](control_print_dof_line.md).

## Example

```
control_print_dof_line           40  -separate_index
control_print_dof_line_coordinates 40  0.5 0.5 0.5 0.6
control_print_dof_line_n         40  2
control_print_dof_line_move      40  -yes
control_timestep                 40  0.1 0.2
```
