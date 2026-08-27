# control_print_dof_smooth_n

## Description

`control_print_dof_smooth_n` sets the NUMBER of smoothing passes of
[`control_print_dof_smooth_dof`](control_print_dof_smooth_dof.md) with
the same index (manual Professional 6.272). If it is not specified, the
smoothing is done **10 times**.

More passes give increasingly smooth (more averaged) results; repeated
passes converge towards the mean value of the smoothed dof.

## Uso

```
control_print_dof_smooth_dof 0  -velx
control_print_dof_smooth_n   0  3
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `0`       | Index of the control record. Must match `control_print_dof_smooth_dof`. |
| `number_of_smoothings` | Number of smoothing passes (>= 1). |

## Output

No output by itself; it modifies the values written to `dof.<index>` by
[`control_print_dof`](control_print_dof.md).

## Example

With the chain of [`control_print_dof_smooth_dof`](control_print_dof_smooth_dof.md)
(`0, 1, 2, 3, 4`), 3 passes give `1.5, 1.5, 2, 2.5, 2.5`; the default 10
passes give `1.9375, 1.96875, 2, 2.03125, 2.0625` (converging to the
mean value 2).

## Differences with the Professional version

- None. `number_of_smoothings < 1` is an error.
