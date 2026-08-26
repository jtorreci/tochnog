# control_print_history_relative_time

## Description

`control_print_history_relative_time` shifts the time axis of the
history files: the time printed in the history files is not the actual
time but the actual time minus the relative time `tr`. This is
convenient when something happens suddenly after a long time — the long
initial time is removed from the time axis of the history plot, making
it clear.

It works in combination with the `control_print_history` record (same
index). The smoothed history files of
`control_print_history_smooth` share the same relative time axis.

## Uso

```
control_print_history                 20  -node_dof 1 -disx
control_print_history_relative_time   20  0.2
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match the `control_print_history` index. |
| `tr`      | Relative time (a real number, in the time units of the calculation) subtracted from the printed time. |

## Output

The history files (`<dof><index>.his`, `smooth...` files) keep the same
layout; only the first column (time) changes: `time_current - tr`.

## Example

```
control_print_history               20  -node_dof 1 -disx
control_print_history_relative_time 20  0.2
control_timestep                    20  0.1 0.3
```

With a calculation running from t = 0.1 to t = 0.3, the history file
prints the times -0.1, 0.0, 0.1.

## Validation

Test `hreltime` (validation-suite/test-2014): block 20 runs t = 0.1..
0.3 with `tr = 0.2` — the last line of `disx1.his` is 0.1 (0.3-0.2).
Block 21 runs without the record — the last line of `disx2.his` is 0.4.
The 0.3 difference between the two files discriminates the shift.
