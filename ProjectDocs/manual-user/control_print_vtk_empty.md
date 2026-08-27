# control_print_vtk_empty

## Description

`control_print_vtk_empty` controls whether the empty elements are
included in the `.vtk` files of `control_print_vtk` (same `icontrol`
index). `-yes` (default) includes them, `-no` excludes them.

An element is "empty" when the solver marks it with
`element_empty -yes`: this happens automatically for
`materi_diffusion` / `materi_density` models when none of the nodes of
the element has a value above the minimum (`EPS_MATERI_DIFFUSION_MINIMUM`
= 0.5, `EPS_MATERI_DENSITY_MINIMUM` = 1e-9). Empty elements are skipped
by the element loop, so they carry no results.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_print_vtk` record:

```
control_print_vtk            36  -yes
control_print_vtk_empty      36  -no
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `36`      | Index of the control record. Must match a `control_print_vtk` index. |
| switch    | `-yes` (default): empty elements are included in the vtk file. `-no`: empty elements are excluded from `CELLS` and `CELL_TYPES` (the counts stay consistent). |

## Example

A model with two elements where the second one is empty (density 0 on
all its nodes):

```
control_print_vtk            35  -yes
control_print_vtk            36  -yes
control_print_vtk_empty      36  -no
```

`tn35.vtk` (default) contains both cells; `tn36.vtk` contains only the
non-empty element (`CELLS 1 5`, `CELL_TYPES 1`).
