# control_print_materi_stress_force

## Description

`control_print_materi_stress_force` (manual Professional 6.328) prints
the forces and moments calculated by `post_calcul -materi_stress
-force` to a special purpose ASCII file for external postprocessing.
The file name is `materi_stress_force.<index>`, where `index` is the
RECORD index (e.g. `control_print_materi_stress_force 100 -all`
produces `materi_stress_force.100`, exactly the manual example). The
files themselves contain comments explaining the detailed structure.

The method is `-all` (every result, including the averaged ones) or
`-primary` (only the primarily calculated results, i.e. WITHOUT the
quad9/hex27 averaged results of
[`post_calcul_materi_stress_force_average`](post_calcul_materi_stress_force_average.md)).

## Input syntax

```
control_print_materi_stress_force <index> -all | -primary
```

## Output

One line per node:

```
node norx_sig nory_sig nors_sig shex_sig shey_sig shes_sig momx_sig momy_sig moms_sig   (2D)
node norx_sig nory_sig norz_sig nors_sig shex_sig shey_sig shez_sig shes_sig mom1x_sig mom1y_sig mom1z_sig mom1s_sig mom2x_sig mom2y_sig mom2z_sig mom2s_sig   (3D)
```

The `x/y/z` components are GLOBAL PLOT components (the vector is drawn
in the structure thickness direction; they are NOT the physical force
components); the `s` component is the physical vector size (the square
root of the sum of the squares), usable for design.

## Example

```
post_calcul_materi_stress_force_element_group 0
post_calcul -materi_stress -force
control_print_materi_stress_force 100 -all
```

## Notes

- The numerical integration is NOT implemented yet (lot 1 =
  infrastructure). Until lots 2/3 land the printed values are 0; the
  file structure (header comments, one line per node, `-all`/`-primary`
  methods) is functional.
- Without a `post_calcul -materi_stress -force` block no file is
  written (documented decision).
- The file is appended at every print call (one block per call).
