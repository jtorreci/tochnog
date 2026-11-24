# post_calcul_materi_stress_force_direction_exclude_epsilon

## Description

`post_calcul_materi_stress_force_direction_exclude_epsilon` (manual
Professional 6.910) is the tolerance `eps` of the exclusion test of
[`post_calcul_materi_stress_force_direction_exclude`](post_calcul_materi_stress_force_direction_exclude.md):
an element side is excluded when `|n . dir| > 1 - eps`. A small `eps`
excludes only very precise normals; a large `eps` also excludes
imprecise ones. Default `1.e-8`.

## Input syntax

```
post_calcul_materi_stress_force_direction_exclude_epsilon eps
```

## Notes

- Single value, 3D concept (in 2D the direction records are ignored).
- Consumed by the numerical integration (lots 2/3); the default
  `1.e-8` applies when the record is absent.
