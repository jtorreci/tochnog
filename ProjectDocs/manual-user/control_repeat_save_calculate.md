# control_repeat_save_calculate

## Description

`control_repeat_save_calculate` performs a statistical analysis of the
data collected by [`control_repeat_save`](control_repeat_save.md) once
the `control_repeat` loop has finished. For every saved data item it
computes the average value (mean) and the variance, and stores them in
the record `repeat_calculate_result`.

Use it at the end of a Monte-Carlo style repeat loop to obtain the mean
and the variance of each sampled quantity directly in Tochnog, without
external post-processing.

## Uso

Place the record at the SAME `icontrol` index as the `control_repeat`
and `control_repeat_save` records:

```
control_timestep             10  1. 100.
control_print                20  -node_dof 1 -sigxx
control_repeat               30  80 10
control_repeat_save          30  -post_point_dof 1 -disy
control_repeat_save_calculate 30  -yes
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `30`      | Index of the control record. Must match the `control_repeat` index. |
| `switch`  | `-yes` to perform the statistical analysis when the repeat completes. |

## Output

The statistics are stored in `repeat_calculate_result` (DOUBLE), one
index per saved data item, with two values each:

- `repeat_calculate_result[j][0]` = average of data item `j`
- `repeat_calculate_result[j][1]` = variance of data item `j`
  (population variance, sum of squared deviations divided by the
  number of repeats)

## Example

```
control_timestep             10  0.05 0.05
control_repeat               20  3 10
control_repeat_save          20  -post_point_dof 1 -veliy
control_repeat_save_calculate 20  -yes

target_item  1  -repeat_calculate_result 0 0
target_value 1  -0.001 1.e-9
target_item  2  -repeat_calculate_result 0 1
target_value 2  1.6666667e-7 1.e-9
end_data
```

The three saved values `-0.0005, -0.001, -0.0015` (arithmetic
sequence, see the `mrepeat_save` test) give average `-0.001` and
variance `d^2*(N^2-1)/12 = 1.6666667e-7`.

## Validation

Test `mrepeat_save` (validation-suite/test-2014): exact analytic mean
and variance over a known arithmetic sequence. A/B: without this
record, `repeat_calculate_result` does not exist.
