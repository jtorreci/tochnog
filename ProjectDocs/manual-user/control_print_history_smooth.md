# control_print_history_smooth

## Description

`control_print_history_smooth` smooths the data values printed by
`control_print_history` (same `icontrol` index) by computing a moving
average of the last N values for each data value. The smoothed results are
written to separate history files whose names start with `smooth`.

Use this to filter out high-frequency noise from time series such as
stresses, displacements or void ratio when plotting their evolution.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_print_history` record:

```
control_print_history            20  -node_dof 1 -hisv0 -node_dof 1 -sigyy
control_print_history_smooth     20  3
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match the `control_print_history` index. |
| `3`       | Window size N. One integer per data value of `control_print_history`, or a single integer that applies to all data values. |

With `3`, each data value is replaced by the average of the last 3 values
at that time step (fewer values are averaged near the start of the
calculation, where not enough data is available yet).

## Output

One history file per data value, with the same naming scheme as
`control_print_history` but prefixed with `smooth`:

- `smooth<number><index>.his` for plain numbers
- `smooth<dofname><index>.his` for dof labels (e.g. `-hisv0` -> `smoothhisv01.his`)

Each line contains `time value` (the averaged value).

## Example

```
control_timestep                20  0.001 0.04
control_print_history           20  -node_dof 1 -hisv0 -node_dof 1 -sigyy
control_print_history_smooth    20  3
```

Produces `hisv01.his` / `sigyy1.his` (raw) and
`smoothhisv01.his` / `smoothsigyy1.his` (smoothed, moving average of 3).
