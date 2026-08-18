# slide_axisymmetric

## Description

`slide_axisymmetric` marks the slide friction defined by `slide_geometry`
as **axisymmetric**. In an axisymmetric analysis the slide friction acts on
the whole ring of circumference `2*pi*r` around the symmetry axis, where
`r` is the radial distance of the node to the axis (the x-coordinate). The
friction force computed by the slide law is multiplied by `2*pi*r`.

Set the switch to `-yes` to activate the axisymmetric scaling.

Requires `slide_geometry` and `slide_friction` (the Coulomb friction law).

## Uso

Place it in the data part, together with the slide records:

```
slide_geometry     0  -geometry_line 1
slide_friction     0  0.5
slide_axisymmetric 0  -yes
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `0`       | Index of the slide record (must match `slide_geometry`). |
| switch    | `-yes`: axisymmetric slide friction (scaled by `2*pi*r`). `-no` (default): no scaling. |

## Related

- `slide_geometry` — the geometry on which the sliding occurs.
- `slide_friction` — Coulomb friction parameters.
- `group_axisymmetric` — the element-group flag that marks the analysis as
  axisymmetric.

## Estado de implementación

- **Implementado**: `slide_axisymmetric` (scale the slide friction by
  `2*pi*r`, r = node x-coordinate). Validated with `slide_axi` (axisymmetric
  quad4 with slide friction active, the model converges and the slide node is
  braked).
