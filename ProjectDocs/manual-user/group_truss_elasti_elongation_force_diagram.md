# group_truss_elasti_elongation_force_diagram

## Description

`group_truss_elasti_elongation_force_diagram index l_0 F_0 l_1 F_1 ...`
(manual 6.773, index = element group): truss force F as a function of
l = elongation / initial truss length.

## Effect on the calculation

- The force is read from the diagram at the current total l (linear
  interpolation) and replaces the elastic force. It is nonlinear elastic:
  on unloading the diagram is followed back.
- Outside the diagram the force is 0.
- `group_truss_elasti_young` and `group_truss_initial_force` are ignored
  when the diagram is given (the manual says they cannot be combined; the
  Professional does not stop).
- `group_truss_plasti` and `group_truss_rope` still limit the force.
- It works with `materi_displacement` and with
  `materi_velocity_integrated`.
- An odd number of values, or fewer than two pairs, stops GNU with an
  error (the Professional gives force 0 without a message).

Measured on the Professional 02-08-2026, truss length 2, diagram
`-0.01 -20  0 0  0.01 10  0.02 15`:

| elongation | l | force |
|---|---|---|
| 0.01 | 0.005 | 5 |
| 0.028 | 0.014 | 12 |
| 0.04 | 0.02 | 15 |
| 0.05 | 0.025 | 0 |
| -0.01 | -0.005 | -10 |
| -0.04 | -0.02 | 0 |

## Example

```
group_type 0 -truss
group_truss_area 0 1.
group_truss_elasti_elongation_force_diagram 0  -0.01 -20  0 0  0.01 10  0.02 15
```

Self-checking inputs in `validation-suite/test-2026`:
`truss_force_diagram1` .. `truss_force_diagram5`,
`truss_force_diagram_error1`.
