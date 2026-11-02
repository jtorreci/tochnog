# control_print_vtk_coord

## Description

`control_print_vtk_coord` controls whether the node coordinates are
written as a data field to the `.vtk` files of `control_print_vtk` (same
`icontrol` index). `-yes` (default) adds the point field `VECTORS coord`
(the stored node coordinates) at the start of `POINT_DATA`; `-no` drops
only that field.

The `POINTS` block of the mesh is always written, whatever the switch: a
vtk file without points could not be opened by a post-processor. With
`control_print_vtk_node_type -node_deformed_mesh` the `POINTS` hold the
deformed coordinates, while `coord` keeps the stored (undeformed) ones.

The global [print_vtk_coord](print_vtk_coord.md) gives its default for every index.

## Usage

Place it in the data part, with the same `icontrol` index as the
`control_print_vtk` record:

```
control_print_vtk            30  -yes
control_print_vtk_coord      30  -no
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `30`      | Index of the control record. Must match a `control_print_vtk` index. |
| switch    | `-yes` (default): the `VECTORS coord` point field is written. `-no`: that field is omitted; `POINTS` stays. |

## Example

```
control_print_vtk            30  -yes
control_print_vtk_coord      30  -no
```

Produces a `.vtk` file with the full mesh (`POINTS`, `CELLS`,
`CELL_TYPES`) and the field data, without the `coord` point field.

## Differences with the Professional version

None: same meaning as the Professional 02-08-2026 (BACKLOG C98, resolved
2026-11-02). Before F5 the GNU `-no` omitted the `POINTS` block and the
GNU had no `coord` field.
