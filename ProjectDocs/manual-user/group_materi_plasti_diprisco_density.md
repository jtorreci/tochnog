# group_materi_plasti_diprisco_density

## Description

di Prisco plasticity with varying density (manual Professional 6.696 and
the theory section "di Prisco plasticity model with varying density"). It is
the `group_materi_plasti_diprisco` model, but with two parameter sets, one
for loose soil and one for dense soil. The ten di Prisco parameters actually
used are interpolated linearly between the two sets with the current void
ratio `e`:

```
w = (el - e) / (el - ed)        clamped to [0, 1]
p = (1 - w) p_loose + w p_dense  (for each of the ten parameters)
```

So `e >= el` uses the loose set, `e <= ed` the dense set. The current void
ratio is history variable 12 (`dipriscohis11`). Give it its initial value
(for example with `control_reset_dof -dipriscohis11`); it then evolves with
the total volumetric strain increment of every step:

```
e_new = e + (1 + e) tr(d_eps)
```

The interpolation uses the void ratio at the start of the step.

## Syntax

```
group_materi_plasti_diprisco_density index gl bfl bpl cpl tpl tcl tel xcl xel bf0l
                                           gd bfd bpd cpd tpd tcd ted xcd xed bf0d el ed
```

| Parameter | Meaning |
|---|---|
| first ten | loose soil di Prisco data, as in `group_materi_plasti_diprisco` (gamma, beta_hat_f, bp, cp, tp, theta_hat_c, theta_hat_e, xi_c, xi_e, beta_f0) |
| next ten | dense soil di Prisco data, same order |
| `el` | void ratio of the loose soil |
| `ed` | void ratio of the dense soil (`el` should be larger than `ed`) |

Requires `materi_stress`, `materi_strain_plasti` and
`materi_plasti_diprisco_history 12` (the eleven variables of the di Prisco
model plus the void ratio).

Errors:

- `Error in group_materi_plasti_diprisco_density. The el should be larger than the ed.`
- `Error: group_materi_plasti_diprisco_density needs materi_plasti_diprisco_history 12.`

## Example

```
materi_plasti_diprisco_history 12
...
group_materi_plasti_diprisco_density 0 3.62 0.5 0.0049 18. 10. 0.253 0.1184 -0.2585 -0.1178 1.1
  4.5 0.6 0.004 25. 12. 0.3 0.15 -0.3 -0.15 1.2 0.9 0.6
control_reset_dof            5 -dipriscohis11
control_reset_value_constant 5 0.75
```

Tests: `validation-suite/test-2026/group_materi_plasti_diprisco_density1`-`4`,
`_error1`, `_error2`.

## Differences with the Professional version

- The Professional 02-08-2026 uses only the first (loose) parameter set:
  the dense set and the value of history variable 12 change nothing, and
  the variable is stored unchanged. This contradicts its manual ("the
  actually used data will be interpolated between the loose and dense data
  using the current density"). The GNU follows the manual (DIVERGENCES
  DIPRISCO-DENSITY, class D, owner decision 2026-11-07).
- The manual calls variable 12 the "relative density (for example 20 or
  40)", but the Professional requires `el > ed`, so `el` and `ed` are void
  ratios; the GNU reads variable 12 as the current void ratio.
