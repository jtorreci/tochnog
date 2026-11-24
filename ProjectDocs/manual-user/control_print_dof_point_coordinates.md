# control_print_dof_point_coordinates

## Description

`control_print_dof_point_coordinates` gives the coordinates of the point
used by [`control_print_dof_point`](control_print_dof_point.md) with the
same index (manual Professional 6.282).

In 1D only the x-coordinate is needed, in 2D `x y`, in 3D `x y z`.

## Uso

```
control_print_dof_point_coordinates 44  0.5 0.5
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `44`      | Index of the control record. Must match the `control_print_dof_point` record. |
| `x y z`   | Coordinates of the point (`ndim` values). |

## Output

No output by itself; it defines the point printed by
[`control_print_dof_point`](control_print_dof_point.md).

## Example

```
control_print_dof_point            44  -separate_index
control_print_dof_point_coordinates 44  0.5 0.5
```
