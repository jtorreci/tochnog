# control_print_vtk_coord

## Description

`control_print_vtk_coord` controls whether the node coordinates are
written to the `.vtk` files of `control_print_vtk` (same `icontrol`
index). `-yes` (default) writes the coordinates, `-no` does not.

In the current VTK writer the coordinates appear only in the `POINTS`
block of the mesh (there are no coordinate data fields like the gmsh
`node_0`/`node_1`/`node_2` scalars), so `-no` omits the whole `POINTS`
block while keeping `CELLS`/`CELL_TYPES` and the `POINT_DATA` fields.

> Note: a VTK file without `POINTS` is not a valid dataset for a
> post-processor (the cells reference point indices that do not exist).
> The switch is meant for smaller debug dumps / for extracting only the
> field data; keep the default `-yes` for visualisation.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_print_vtk` record:

```
control_print_vtk            30  -yes
control_print_vtk_coord      30  -no
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `30`      | Index of the control record. Must match a `control_print_vtk` index. |
| switch    | `-yes` (default): node coordinates (POINTS block) are written. `-no`: the POINTS block is omitted. |

## Example

```
control_print_vtk            30  -yes
control_print_vtk_coord      30  -no
```

Produces a `.vtk` file with the mesh topology and the field data but
without the coordinate block. Remove the `control_print_vtk_coord`
record (or use `-yes`) to restore the default output.
