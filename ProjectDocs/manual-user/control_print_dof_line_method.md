# control_print_dof_line_method

## Description

`control_print_dof_line_method` selects which node coordinates are used
for the interpolation AND written to the files of
[`control_print_dof_line`](control_print_dof_line.md) (manual
Professional 6.277):

- `-node`: the values of the `node` record (the CURRENT coordinates).
- `-node_start_refined`: the values of the `node_start_refined` record
  (the REFERENCE coordinates).

In an updated Lagrange formulation where the mesh nodes follow the
material the two differ; in a geometrically linear analysis they do not.

Default: `-node_start_refined` (falls back to the stored node
coordinates when no `node_start_refined` record exists).

## Uso

```
control_print_dof_line_method 38  -node_start_refined
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `38`      | Index of the control record. Must match the `control_print_dof_line` record. |
| `node_type` | `-node` or `-node_start_refined` (default). |

## Output

No output by itself; it selects the coordinate frame of the line
geometry. The printed coordinates of each point are the line positions
in that frame; the values are interpolated with the shape functions of
the element in that frame.

## Example

```
control_print_dof_line           38  -separate_index
control_print_dof_line_coordinates 38  0.5 0. 0.5 0.995
control_print_dof_line_n         38  2
control_print_dof_line_method    38  -node_start_refined
```
