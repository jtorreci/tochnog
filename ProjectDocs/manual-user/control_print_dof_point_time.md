# control_print_dof_point_time

## Description

`control_print_dof_point_time` writes the `time_current` as the FIRST
line of each file of [`control_print_dof_point`](control_print_dof_point.md)
with the same index, in gnuplot comment format (manual Professional
6.283): `# time <time_current>`.

The comment is written once per call (each time the file is appended).

## Uso

```
control_print_dof_point_time 43  -yes
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `43`      | Index of the control record. Must match the `control_print_dof_point` record. |
| switch    | `-yes`: write the time comment. `-no` (default): do not. |

## Output

No output by itself; it prepends `# time <time_current>` to the data
lines of [`control_print_dof_point`](control_print_dof_point.md).

## Example

```
control_print_dof_point            43  -separate_index
control_print_dof_point_coordinates 43  0.5 0.5
control_print_dof_point_time       43  -yes
control_timestep                   43  0.1 0.1
```

Produces `disy.43`:

```
# time 0.1
0.5 0.5 ...
```
