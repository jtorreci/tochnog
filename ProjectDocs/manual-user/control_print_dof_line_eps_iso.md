# control_print_dof_line_eps_iso

## Description

`control_print_dof_line_eps_iso` sets the tolerance with which a point of
the line of [`control_print_dof_line`](control_print_dof_line.md) is
accepted as part of an element (manual Professional 6.276).

The default value is `1.e-3`. Increase it when the mesh is not exactly
adjusted to the line (e.g. a line slightly outside the mesh boundary).

## Uso

```
control_print_dof_line_eps_iso 37  1.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `37`      | Index of the control record. Must match the `control_print_dof_line` record. |
| `eps_iso` | Tolerance (positive). Default `1.e-3`. |

## Output

No output by itself; it controls whether each line point is accepted in
the point-in-element test of
[`control_print_dof_line`](control_print_dof_line.md). Points that are
not accepted are omitted from the printed files.

## Example

```
control_print_dof_line           37  -separate_index
control_print_dof_line_coordinates 37  0.5 0.5 1.02 0.5
control_print_dof_line_n         37  2
control_print_dof_line_eps_iso   37  1.0
```

The point `(1.02, 0.5)` is 0.02 outside a mesh spanning `x in [0,1]`;
with the default tolerance it is omitted, with `eps_iso = 1.0` it is
accepted (extrapolated values).
