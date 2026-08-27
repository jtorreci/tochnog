# post_calcul_materi_stress_force_plot_switch

## Description

`post_calcul_materi_stress_force_plot_switch` (manual Professional
6.916) inverts the direction in which the vectors are drawn: set the
switch to `-yes` for the vectors you want drawn the other way. In 2D
you need one switch for the normal force, the shear force and the
moment (3 switches); in 3D one switch for the normal force, the shear
force and the two moments (4 switches).

## Input syntax

```
post_calcul_materi_stress_force_plot_switch switch_0 switch_1 ...
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `switch_i` | `-yes` inverts the vector direction, `-no` (default) keeps it. One switch per vector item: 3 in 2D (nor, she, mom), 4 in 3D (nor, she, mom1, mom2). |

## Notes

- The number of switches must equal 3 (2D) or 4 (3D); otherwise an
  error is raised.
- Consumed by the numerical integration (lots 2/3): the record is
  validated and stored.
