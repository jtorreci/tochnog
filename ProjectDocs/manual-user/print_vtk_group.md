# print_vtk_group

## Description

`print_vtk_group switch` (manual Professional 6.1001): with `-yes`
(default) every vtk file (`control_print_vtk`, `print_vtk_calculation`)
holds the element group of each cell as the cell field `element_group`
(`CELL_DATA` block, group 0 for an element without `element_group`
record), so that paraview can filter the groups with a Threshold on
`element_group`. With `-no` the field is not written.

## Usage

```
print_vtk_group -no
```

| Parameter | Meaning |
|---|---|
| switch | `-yes` (default): write the `element_group` cell field. `-no`: do not. |

The field follows the cells written (`control_print_vtk_empty -no` drops
the empty elements from both). Tests `validation-suite/test-2026/print_vtk_group1/2`.

## Differences with the Professional version

- The Professional 02-08-2026 ignores `-no` and always writes
  `element_group`; the GNU follows the manual (DIVERGENCES `VTK-GROUP-NO`,
  class D).
- Before F5 the GNU vtk files had no cell field; the field is new output in
  every GNU vtk file (manual default, = Professional).
