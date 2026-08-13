# control_reset_dof

## Description

`control_reset_dof` resets the value of a nodal dof to a target value
that is either constant or depends on another dof. This is useful, for
example, to make the void ratio (`hisv0`) of a hypoplastic law depend on
the vertical stress (`sigyy`).

The target value is defined by one of:

- `control_reset_value_constant` — a constant value.
- `control_reset_value_dof` + `control_reset_value_dof_diagram` — a value
  taken from another dof via a table (linear interpolation).

The method is set by `control_reset_value_method`:
`-use` (default) sets the dof to the value, `-add` adds it,
`-multiply` multiplies the current value.

## Uso

Place it in the data part:

```
control_reset_dof                20  -hisv0
control_reset_value_constant     20  0.55
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. |
| `-hisv0`  | The nodal dof to reset. |

## Related records

- `control_reset_value_constant index value` — constant target value.
- `control_reset_value_dof index dof` — target dof whose value drives the
  reset.
- `control_reset_value_dof_diagram index z_0 value_0 z_1 value_1 ...` —
  table (dof value → target value), linear interpolation.
- `control_reset_value_method index method` — `-use` (set), `-add`,
  `-multiply`.

## Example

Make the void ratio (`hisv0`) of the hypo law depend on the vertical
stress (`sigyy`):

```
control_reset_dof                10  -hisv0
control_reset_value_dof          10  -sigyy
control_reset_value_dof_diagram  10  -1.e2 0.4  -2.e2 0.38  -1.e3 0.3
control_reset_value_method       10  -use
```

For vertical stress -100 the void ratio is reset to 0.4, for -200 to
0.38, etc. (linear interpolation between points).
