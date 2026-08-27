# post_calcul_materi_stress_force_average

## Description

`post_calcul_materi_stress_force_average` (manual Professional 6.908)
is only available for quad9 and hex27 elements. When forces and moments
are primarily calculated in two opposing end faces of the element, with
`-yes` the forces and moments of the nodes in the plane between the two
end faces are set to the averaged values of the two opposing end faces.
Default `-yes`.

## Input syntax

```
post_calcul_materi_stress_force_average -yes | -no
```

## Notes

- Single switch, consumed by the numerical integration (lots 2/3). The
  `-primary` method of
  [`control_print_materi_stress_force`](control_print_materi_stress_force.md)
  prints only the non-averaged (primary) results.
- Until the averaging exists (lots 2/3), every node is primary, so
  `-all` and `-primary` print the same lines.
