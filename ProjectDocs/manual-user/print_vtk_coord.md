# print_vtk_coord

## Description

`print_vtk_coord switch` (manual Professional 6.1000): default of
[control_print_vtk_coord](control_print_vtk_coord.md) for every vtk file. A
`control_print_vtk_coord` record for a control index overrules it for that
index. In the GNU the node coordinates are the `POINTS` block, so `-no`
omits it (same meaning as `control_print_vtk_coord -no`). Default `-yes`.

## Usage

```
print_vtk_coord -no
control_print_vtk_coord 3 -yes   # index 3 keeps its POINTS
```

| Parameter | Meaning |
|---|---|
| switch | `-yes` (default): write the coordinates. `-no`: omit them. |

Tests `validation-suite/test-2026/print_vtk_coord1/2`.

## Differences with the Professional version

The Professional 02-08-2026 always writes `POINTS` and also a point field
`coord` (the node coordinates); `print_vtk_coord -no` and
`control_print_vtk_coord -no` drop that field. The GNU has no `coord`
field and drops `POINTS` (BACKLOG C98, DIVERGENCES `C98`, class B, open:
owner decision on the meaning of the switch).
