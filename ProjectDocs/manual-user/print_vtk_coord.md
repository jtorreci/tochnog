# print_vtk_coord

## Description

`print_vtk_coord switch` (manual Professional 6.1000): default of
[control_print_vtk_coord](control_print_vtk_coord.md) for every vtk file. A
`control_print_vtk_coord` record for a control index overrules it for that
index. `-yes` (default) adds the point field `VECTORS coord` (the stored
node coordinates) to every vtk file; `-no` drops only that field, the
`POINTS` block is always written.

## Usage

```
print_vtk_coord -no
control_print_vtk_coord 3 -yes   # index 3 keeps its coord field
```

| Parameter | Meaning |
|---|---|
| switch | `-yes` (default): write the `coord` point field. `-no`: omit it. |

Tests `validation-suite/test-2026/print_vtk_coord1/2`.

## Differences with the Professional version

None: same meaning as the Professional 02-08-2026 (BACKLOG C98, resolved
2026-11-02 by owner decision).
