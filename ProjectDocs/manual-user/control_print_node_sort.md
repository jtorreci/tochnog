# control_print_node_sort

## Description

`control_print_node_sort` sorts the lines printed by
[`control_print_node`](control_print_node.md) with the same index
(manual Professional 6.334), in ASCENDING order of the sort key:

- with [`control_print_node_angular`](control_print_node_angular.md):
  `-angle` (sort by the printed angle);
- otherwise: `-x`, `-y` (2D/3D) or `-z` (3D) — sort by that coordinate.

Without the record the lines are written in node order. The manual's
mention of "`control_print_node_method`" is a typo for
`control_print_node_angular`; it is not a keyword.

## Uso

```
control_print_node            57  -node_dof -vely
control_print_node_sort       57  -y
control_timestep             57  1.0 1.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `57`      | Index of the control record. Must match `control_print_node`. |
| `sort_method` | `-angle` (with angular), `-x`, `-y` (2D/3D), `-z` (3D). |

## Output

Modifies the files of [`control_print_node`](control_print_node.md):
the lines are written sorted ascending by the sort key (stable order
for equal keys).

## Example

A square whose node NUMBERS are not ordered by `y` (node 1 at `(0,1)`,
node 2 at `(0,0)`, ...): without sort the first line is `0 1 -0.01`;
with `-y` the first line is `0 0 0` and the last `1 1 -0.01`.

## Differences with the Professional version

- None. `-angle` without angular, `-y` in 1D and `-z` in 2D are errors.
