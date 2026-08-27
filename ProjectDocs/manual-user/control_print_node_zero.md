# control_print_node_zero

## Description

`control_print_node_zero` suppresses or activates the printing of
results with value zero by
[`control_print_node`](control_print_node.md) with the same index
(manual Professional 6.335). The default (record not given) is `-yes`:
zero valued results ARE printed. With `-no` the lines whose printed
value is zero are omitted (per file, i.e. per selected part).

## Uso

```
control_print_node      50  -node_dof -vely
control_print_node_zero 50  -no
control_timestep       50  1.0 1.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `50`      | Index of the control record. Must match `control_print_node`. |
| switch    | `-yes` (default): print zero valued results. `-no`: suppress them. |

## Output

Modifies the files of [`control_print_node`](control_print_node.md):
lines with a zero value are omitted for the parts where the value is
zero.

## Example

A model with `vely = 0` at the bottom and `-0.01` at the top: with
`-no` the `vely` file has only the 2 lines of `-0.01`; a `velx` file
(velx = 0 everywhere) becomes empty.

## Differences with the Professional version

- The zero comparison is EXACT (`value == 0.0`).
