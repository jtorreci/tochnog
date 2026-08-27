# control_print_node_geometry

## Description

`control_print_node_geometry` restricts
[`control_print_node`](control_print_node.md) with the same index to
print ONLY the nodes located on the given geometry (manual Professional
6.333). A node is considered on the geometry with the same test used by
the `print_filter` records (projection of the node onto the geometry
entity within its tolerance).

## Uso

```
geometry_line 0  0. 0. 1. 0. 1.e-5

control_print_node        54  -node_dof -vely
control_print_node_geometry 54  -geometry_line 0
control_timestep         54  1.0 1.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `54`      | Index of the control record. Must match `control_print_node`. |
| `geometry_item_name` | Geometry item, e.g. `-geometry_line`, `-geometry_point`, `-geometry_circle`, ... |
| `geometry_item_index` | Index of the geometry item. |

## Output

Modifies the files of [`control_print_node`](control_print_node.md):
nodes outside the geometry are omitted.

## Example

With `geometry_line 0` on the bottom edge `(0,0)-(1,0)` of a unit
square, only the 2 bottom nodes are printed (2 lines instead of 4).

## Differences with the Professional version

- The node is projected with the `-node_start_refined` frame when such
  records exist, falling back to the stored node coordinates (same
  fallback as `control_print_dof_line`).
