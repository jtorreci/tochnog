# control_repeat_save

## Description

`control_repeat_save` selects data items whose current value must be
saved every time a `control_repeat` repeat is performed. Each repeat
writes one record `repeat_save_result` (index 0 for the first repeat,
index 1 for the second, and so on) holding the values of all selected
items at the moment of the jump back to the repeat start.

Use it together with `control_repeat` to collect one scalar (or a few
scalars) per realization of a Monte-Carlo style analysis — e.g. the
final displacement of a point, a stress component, or the current
time — and then post-process the collected values (statistics can be
computed in Tochnog itself with
[`control_repeat_save_calculate`](control_repeat_save_calculate.md),
or externally with the `.dbs` dump / SQLite / CSV output).

## Uso

Place the record at the SAME `icontrol` index as the `control_repeat`
record, in the data part:

```
control_timestep             10  1. 100.
control_print                20  -node_dof 1 -sigxx
control_repeat               30  80 10
control_repeat_save          30  -post_point_dof 1 -disy
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `30`      | Index of the control record. Must match the `control_repeat` index. |
| `data_item_name_0 data_item_index_0 data_item_number_0 ...` | List of triplets, one per data item to save: the record name (e.g. `-post_point_dof`, `-time_current`), its index, and the value number within it (a plain number, or a negative dof label such as `-disy` / `-sigxx`). |

## Output

For every repeat the current values are stored in the record
`repeat_save_result` (DOUBLE, one index per repeat):

- `repeat_save_result[0]` = values at the moment of the first jump back
  (first repeat), `repeat_save_result[1]` = values of the second, etc.
- Within one index the values appear in the order of the triplets.

The values can be read with `target_item` (as in the validation test)
or via the `.dbs` database dump.

## Example

```
control_timestep             10  0.05 0.05
control_repeat               20  3 10
control_repeat_save          20  -post_point_dof 1 -veliy

target_item  1  -repeat_save_result 0 0
target_value 1  -0.0005 1.e-10
end_data
```

The column is compressed one step per repeat (dt = 0.05, vely = -0.01),
so the integrated vertical displacement saved at each jump is
exactly `-0.0005`, `-0.001`, `-0.0015` (see the `mrepeat_save` test).

## Validation

Test `mrepeat_save` (validation-suite/test-2014): analytic sequence and
exact targets on `repeat_save_result` and on the statistics computed by
`control_repeat_save_calculate`. A/B: without `control_repeat_save_calculate`
the record `repeat_calculate_result` does not exist (a target on it
fails).
