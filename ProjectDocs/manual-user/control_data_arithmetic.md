# control_data_arithmetic

## Description

Change a data item arithmetically during the calculation. With
`data_item_name data_item_index data_item_number` you select which data
item to change; it is changed with the value `val` specified in the
corresponding `control_data_arithmetic_double` record (same index).
With `operat` you select how: `-plus`, `-minus`, `-multiply` or
`-divide` (manual Professional 6.115/6.116).

Instead of a specific index you can also specify a range
`-ra ... -ra`. In case you specify `-all` for `data_item_number`, the
value is used for ALL numbers of the record.

The record fires at the end of every timestep of the `control_timestep`
block with the same index (the operation is applied once per step).

## Usage

```
control_data_arithmetic <index> <data_item_name> <index|-ra ... -ra> <number|-all> <-plus|-minus|-multiply|-divide>
control_data_arithmetic_double <index> <val>
```

## Example

```
control_timestep 0  0.1 0.1

control_data_arithmetic         1  -group_materi_elasti_young 0 0 -multiply
control_data_arithmetic_double  1  2.0

control_timestep 1  0.1 0.1
```

During timestep block 1 the Young's modulus of group 0 is multiplied by
2.0 at each step (1000 -> 2000 after one step). With number `-all` the
operation applies to every value of the record. Division by zero is an
input error; integer-typed records cannot be changed (use
`change_dataitem` or `control_data_put` for those).

See also `change_dataitem` (time-table driven changes) and
`control_data_copy`.
