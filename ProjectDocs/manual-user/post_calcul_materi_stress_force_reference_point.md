# post_calcul_materi_stress_force_reference_point

## Description

`post_calcul_materi_stress_force_reference_point` (manual Professional
6.914) gives the approximate middle point of the tunnel (or a point at
a large perpendicular distance from a sheet pile) so that Tochnog can
orient the calculated forces and moments consistently outwards/inwards
in the structure thickness direction. One point per element group of
[`post_calcul_materi_stress_force_element_group`](post_calcul_materi_stress_force_element_group.md)
is required: 3 values (x y z) per point in 3D, 2 values (x y) in 2D.

In 3D the record is REQUIRED (error when absent). In 2D, when absent,
a warning is issued and the documented default reference point (0,0) is
used.

## Input syntax

```
post_calcul_materi_stress_force_reference_point x_0 y_0 [z_0] x_1 y_1 [z_1] ...
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `x_i y_i [z_i]` | Reference point of the i-th element group (`ndim` values per group). |

## Notes

- The total length must be `ngroups * ndim` (one point per group);
  otherwise an error is raised.
- Consumed by the numerical integration (lots 2/3), not yet
  implemented: the record is validated and stored.
