# force_element_edge_multi_linear_factor_x

## Description

`force_element_edge_multi_linear_factor_x` defines a **multilinear
multiplication factor** in the x-direction that multiplies the force of a
`force_element_edge` record with the same index. This lets you apply
x-coordinate dependent distributed edge forces.

Outside the specified x-range the factor is taken as 0.

## Uso

Place it in the data part, with the same `icontrol` index as the
`force_element_edge` record:

```
force_element_edge          4  0.0   -5.00
force_element_edge_geometry 4  -geometry_line 3
force_element_edge_multi_linear_factor_x 4  0.0 1.0  50.0 1.0  50.0001 0.0  100.0 0.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `4`       | Index of the force record. Must match a `force_element_edge` index. |
| `x_0 factor_0 x_1 factor_1 ...` | Multilinear diagram: pairs of x-coordinate and multiplication factor. Linear interpolation between consecutive points; factor 0 outside the range. |

## Example

Load the top edge of a beam (y = 2.5) with a distributed vertical force
that only acts on the left half (x in [0, 50]):

```
force_element_edge          4  0.0   -5.00
force_element_edge_geometry 4  -geometry_line 3
force_element_edge_multi_linear_factor_x 4  0.0 1.0  50.0 1.0  50.0001 0.0  100.0 0.0
```

The force is 5 N/mm for x in [0, 50] and 0 for x > 50.
