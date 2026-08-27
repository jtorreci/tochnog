# control_print_node_angular

## Description

`control_print_node_angular` makes
[`control_print_node`](control_print_node.md) with the same index write
an ANGLE instead of the coordinates (manual Professional 6.331). The
angle is in **degrees**, measured from the positive global axis directed
to another positive global axis:

| switch_x switch_y switch_z | Angle from | Formula |
|----------------------------|------------|---------|
| `-yes -yes -no`            | +x to +y   | atan2(y−ym, x−xm) |
| `-no -yes -yes`            | +y to +z   | atan2(z−zm, y−ym) |
| `-yes -no -yes`            | +x to +z   | atan2(z−zm, x−xm) |

The middle point `(xm, ym, zm)` comes from
[`control_print_node_angular_middle`](control_print_node_angular_middle.md)
(default `0 0 0`).

In 1D this record cannot be used. In 2D only `-yes -yes` is allowed and
`switch_z` must NOT be given.

## Uso

```
control_print_node            20  -node_dof -vely
control_print_node_angular    20  -yes -yes
control_print_node_angular_middle 20  0.5 0.5
control_timestep             20  1.0 1.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match `control_print_node`. |
| switch_x / switch_y / switch_z | `-yes`/`-no` selecting the axis pair (see the table above). |

## Output

Modifies the files of [`control_print_node`](control_print_node.md):
lines `angle <value>` instead of `x y z <value>`.

## Example

2D with middle `(0.5, 0.5)`: the node `(1, 1)` is at
`atan2(0.5, 0.5) = 45` degrees:

```
control_print_node            20  -node_dof -vely
control_print_node_angular    20  -yes -yes
control_print_node_angular_middle 20  0.5 0.5
control_timestep             20  1.0 1.0
```

Produces `vely.20` with lines `angle value`: `-135 0 / -45 0 /
135 -0.01 / 45 -0.01` (nodes `(0,0) (1,0) (0,1) (1,1)`).

## Differences with the Professional version

- None. The angle is in degrees; the 2D/3D switch combinations are
  validated as described in the manual.
