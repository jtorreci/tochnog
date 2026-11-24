# post_calcul_materi_stress_force_direction_exclude

## Description

`post_calcul_materi_stress_force_direction_exclude` (manual
Professional 6.909) tells Tochnog for which element sides forces and
moments should NOT be determined: all element sides whose normal points
in the given direction are neglected. Typically, in a tunnel
calculation you take the tunnel length direction as `dir_x dir_y dir_z`.

In 3D, either this record or
[`post_calcul_materi_stress_force_direction_include`](post_calcul_materi_stress_force_direction_include.md)
is REQUIRED; the two are mutually exclusive (manual 6.913). In 2D the
record is meaningless: a warning is issued and it is ignored.

## Input syntax

```
post_calcul_materi_stress_force_direction_exclude dir_x dir_y dir_z
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `dir_x dir_y dir_z` | Excluded direction (3 values in 3D). Sides with `|n . dir| > 1 - eps` are excluded, with `eps` from [`post_calcul_materi_stress_force_direction_exclude_epsilon`](post_calcul_materi_stress_force_direction_exclude_epsilon.md) (default `1.e-8`). |

## Notes

- The record length must be exactly `ndim` (3 in 3D), otherwise an
  error is raised.
- The direction is consumed by the numerical integration (lots 2/3),
  not yet implemented: the record is validated and stored.
