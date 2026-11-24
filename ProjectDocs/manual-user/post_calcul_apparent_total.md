# post_calcul_apparent_total

## Description

`post_calcul_apparent_total switch` (manual Professional 6.903). With `-yes`
the apparent Young modulus and Poisson ratio of `post_calcul -materi_stress
-young_apparent / -poisson_apparent` are calculated from the total strains
and total stresses. With `-no` (the default) the increments of the last
time step are used.

The values come from the bulk and shear moduli `K = p/eps_v` and
`G = q/(3 eps_q)`: `E = 9KG/(3K+G)`, `nu = (3K-2G)/(2(3K+G))`.
`materi_strain_total` is needed.

## Usage

```
post_calcul -materi_stress -young_apparent -materi_stress -poisson_apparent
post_calcul_apparent_total -yes
```

## Example

`validation-suite/test-2026/post_calcul_apparent_total1.dat` (copy of the
Professional `ground17`): the stresses are reset to -100, then the
volumetric strain -1.e-5 is applied with E 1.e7, nu 0. With the totals:
E 1.333333e7, nu 1/3 (incremental default: 1.e7 and 0).

## Differences with the Professional version

None measured (same values on the example).
