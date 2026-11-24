# control_data_activate

## Description

Set data items to become activated (`-yes`) or de-activated (`-no`)
during the calculation (manual Professional 6.114). The record lists one
or more data item names followed by the switch.

De-activation deletes all records of the listed items: they stop being
used by the solvers from that point on. Typical use: remove a load
family (`bounda_force`), a boundary condition family or gravity at a
given stage of a staged calculation.

GNU difference with Professional: records from the input file are active
by default, so `-yes` is a no-op, and de-activation is DESTRUCTIVE (the
records are deleted; re-activation would require re-putting them with
`control_data_put`).

## Usage

```
control_data_activate <index> <data_item_name_0> <data_item_name_1> ... <-yes|-no>
```

## Example

```
control_timestep 0  0.1 0.1

control_data_activate 1  -bounda_force -no

control_timestep 1  0.1 0.1
```

A column pulled by a `bounda_force` during timestep block 0 loses the
load at block 1 (all `bounda_force` records are deleted): the
quasi-static equilibrium of the unloaded column returns to zero
displacement, while with the force still active it would stay at the
stretched value.
