# post_calcul_absolute

## Description

`post_calcul_absolute switch` (manual Professional 6.902). With `-yes` every
result of `post_calcul` is made positive (absolute value), for example for
a presentation of the results. Default `-no`.

The absolute value is taken last, after `post_calcul_multiply` and
`post_calcul_limit`.

## Usage

```
post_calcul -materi_strain_total -average
post_calcul_absolute -yes
```

## Example

`validation-suite/test-2026/post_calcul_modify1.dat` and
`post_calcul_modify2.dat`: a 1D bar with the total strain -1; the average
-0.333333 becomes positive after the multiply and the limit.

## Differences with the Professional version

None measured: the same values and the same order (multiply, limit,
absolute) as the Professional 02-08-2026.
