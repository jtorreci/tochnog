# control_print_dof_line_time

## Description

`control_print_dof_line_time` writes the `time_current` as the FIRST
line of each file of [`control_print_dof_line`](control_print_dof_line.md)
with the same index, in gnuplot comment format (manual Professional
6.280): `# time <time_current>`.

The comment is written once per call (each time the file is appended).

## Uso

```
control_print_dof_line_time 42  -yes
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `42`      | Index of the control record. Must match the `control_print_dof_line` record. |
| switch    | `-yes`: write the time comment. `-no` (default): do not. |

## Output

No output by itself; it prepends `# time <time_current>` to the data
lines of [`control_print_dof_line`](control_print_dof_line.md).

## Example

```
control_print_dof_line           42  -separate_index
control_print_dof_line_coordinates 42  0.5 0. 0.5 1.
control_print_dof_line_n         42  2
control_print_dof_line_time      42  -yes
control_timestep                 42  0.1 0.1
```

Produces `disy.42`:

```
# time 0.1
0.5 0 ...
0.5 1 ...
```
