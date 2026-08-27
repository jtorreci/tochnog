# control_print_node_angular_middle

## Description

`control_print_node_angular_middle` sets the MIDDLE point of the axes
used by [`control_print_node_angular`](control_print_node_angular.md)
with the same index (manual Professional 6.332). In 2D only `x_middle
y_middle` are given (no `z_middle`); in 3D the three coordinates. The
default, when the record is not given, is `0 0 0`.

The angle follows from e.g. in 2D `tan(angle) = (y−y_middle)/(x−x_middle)`.

In 1D this record cannot be used.

## Uso

```
control_print_node_angular_middle 20  0.5 0.5
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match `control_print_node_angular`. |
| `x_middle y_middle [z_middle]` | Middle point coordinates (2 values in 2D, 3 in 3D). |

## Output

No output by itself; it defines the reference point of the angle written
by [`control_print_node_angular`](control_print_node_angular.md).

## Example

With the middle at the center `(0.5, 0.5)` of a unit square, the
top-right node `(1, 1)` has angle `atan2(0.5, 0.5) = 45` degrees.

## Differences with the Professional version

- None.
