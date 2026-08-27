# post_calcul_materi_stress_force_direction_include

## Description

`post_calcul_materi_stress_force_direction_include` (manual
Professional 6.911) tells Tochnog which element sides produce forces
and moments: all element sides with normals PERPENDICULAR to the given
direction are neglected. Typically, in a sheet pile calculation you
take the sheet pile height direction as `dir_x dir_y dir_z`.

In 3D, either this record or
[`post_calcul_materi_stress_force_direction_exclude`](post_calcul_materi_stress_force_direction_exclude.md)
is REQUIRED; the two are mutually exclusive (manual 6.913). In 2D the
record is meaningless: a warning is issued and it is ignored.

## Input syntax

```
post_calcul_materi_stress_force_direction_include dir_x dir_y dir_z
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `dir_x dir_y dir_z` | Included direction (3 values in 3D). Sides with `|n . dir| < eps` are excluded, with `eps` from [`post_calcul_materi_stress_force_direction_include_epsilon`](post_calcul_materi_stress_force_direction_include_epsilon.md) (default `1.e-8`). |

## Notes

- The record length must be exactly `ndim` (3 in 3D), otherwise an
  error is raised.
- The direction is consumed by the numerical integration (lots 2/3),
  not yet implemented: the record is validated and stored.
