# control_print_beam_force_moment_switch

## Description

Changes the definition of the beam forces and moments printed by
[`control_print_beam_force_moment`](control_print_beam_force_moment.md)
(manual Professional 6.264). Setting the switch to `-yes` multiplies
all 12 components by `-1`, so you can get exactly the sign definition
you want.

## Input syntax

```
control_print_beam_force_moment_switch <index> -yes
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index`   | Index of the control record. Must match the `control_print_beam_force_moment` index. |
| `switch`  | `-yes` inverts the sign of all printed components; `-no` (the default when the record is absent) leaves them unchanged. |

## Example

```
control_print_beam_force_moment 6 -separate_sequential
control_print_beam_force_moment_coordinates 6 0 -0.5 0 0.5
control_print_beam_force_moment_switch 6 -yes
```

## Differences with the Professional version

- `-no` is accepted explicitly (the manual only defines `-yes`); the
  default without the record is no inversion.
