# post_calcul_materi_stress_force_direction_include_epsilon

## Description

`post_calcul_materi_stress_force_direction_include_epsilon` (manual
Professional 6.912) is the tolerance `eps` of the inclusion test of
[`post_calcul_materi_stress_force_direction_include`](post_calcul_materi_stress_force_direction_include.md):
an element side is excluded when `|n . dir| < eps`. A small `eps`
excludes only normals precisely perpendicular; a large `eps` also
excludes imprecise ones. Default `1.e-8`.

## Input syntax

```
post_calcul_materi_stress_force_direction_include_epsilon eps
```

## Notes

- Single value, 3D concept (in 2D the direction records are ignored).
- Consumed by the numerical integration (lots 2/3); the default
  `1.e-8` applies when the record is absent.
