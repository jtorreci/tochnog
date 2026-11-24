# post_calcul_safety_default

## Description

`post_calcul_safety_default eps value` (manual Professional 6.917). For the
hydraulic safety factors of `post_calcul -materi_stress -safety_piping /
-safety_lifting`: when the absolute value of the pressure in the
denominator is below `eps`, the factor is set to `value` (no division by
zero). Without the record eps is very small (1.e-20) and value is 0.

`post_calcul_safety_maximum` does not cap the default value.

## Usage

```
post_calcul -materi_stress -safety_piping -materi_stress -safety_lifting
post_calcul_safety_default 1.e-3 99.
```

## Example

`validation-suite/test-2026/post_calcul_safety_default1.dat` (Professional
`ground14`): above the phreatic level the pressures are 0, so both factors
are 99 there.

## Differences with the Professional version

None measured (same node values).
