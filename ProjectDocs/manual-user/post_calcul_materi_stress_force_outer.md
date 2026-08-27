# post_calcul_materi_stress_force_outer

## Description

`post_calcul_materi_stress_force_outer` (manual Professional 6.915):
with `-yes` the forces and moments are only calculated for the nodes at
the OUTER sides of the elements (the nodes with the furthest distance
relative to the reference point) - gives nicer vector plots. Default
`-no` - gives nicer contour fill plots.

## Input syntax

```
post_calcul_materi_stress_force_outer -yes | -no
```

## Notes

- Single switch, consumed by the numerical integration (lots 2/3).
