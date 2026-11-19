# group_time_fill

## Description

`group_time_fill index birth_empty birth_filled death` (manual 6.770,
index = element group): like `group_time`, but the elements are slowly
filled with material between `birth_empty` and `birth_filled`.

## Effect on the calculation

- The elements are used from `birth_empty` (excluded) to `death`
  (included), as `group_time birth_empty death`.
- Filling factor f = (t - birth_empty) / (birth_filled - birth_empty),
  0 before and 1 after, with t the time at the END of the step.
- The density is multiplied by f (gravity load and mass).
- The total (pore) pressure acting on the skeleton is multiplied by f
  (groundflow pressure unknown and static phreatic level).
- While f < 1 the plasticity data are ignored; from the step that reaches
  `birth_filled` on, plasticity is active.

Measured on the Professional 02-08-2026, two stacked quad4 (E 1000,
nu 0, density 1, gravity -10), upper one `group_time_fill 1 0.5 2.5 100`:

| t | f | disy y=1 | disy y=2 |
|---|---|---|---|
| 0.5 | 0 | -0.005 | 0 |
| 1.0 | 0.25 | -0.0075 | -0.00375 |
| 1.5 | 0.5 | -0.01 | -0.0075 |
| 2.0 | 0.75 | -0.0125 | -0.01125 |
| 2.5 | 1 | -0.015 | -0.015 |

## Example

```
group_time_fill 1  0.5 2.5 100
```

Self-checking inputs in `validation-suite/test-2026`:
`group_time_fill1` (density), `group_time_fill2` (no plasticity while
filling), `group_time_fill3` (pore pressure).
