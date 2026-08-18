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
- `control_reset_value_linear` / `_exponent` / `_power` / `_square_root` /
  `_logarithmic` / `_logarithmic_second` / `_multi_linear` — a value that
  depends on the SPATIAL coordinates of the node (see below).

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

## Distribuciones espaciales (`control_reset_value_*`)

The reset value can depend on the node coordinates (x, y, z). The
coefficients follow the same ordering as the space dimensions (1D: only
`ax...`; 2D: `ax... ay...`; 3D: `ax... ay... az...`).

| Record | Dependency (per dimension) | Coefs per dim |
|--------|---------------------------|---------------|
| `_linear ax ay az` | `ax x + ay y + az z` | 1 |
| `_power ax bx ay by az bz` | `ax x^bx + ay y^by + az z^bz` | 2 |
| `_square_root ax bx cx ...` | `ax sqrt(bx x + cx x^2)` | 3 |
| `_exponent ax bx cx dx ex ...` | `ax e^(bx + cx x dx + ex x)` | 5 |
| `_logarithmic ax bx cx dx ex ...` | `ax ln(bx + cx x dx + ex x)` | 5 |
| `_logarithmic_second ax bx cx dx ex fx gx ...` | `(ax+bx) e^(cx ln(dx (x+ex)/fx)) + gx` | 7 |
| `_multi_linear z0 value0 z1 value1 ...` | multilinear table vs the vertical coordinate (1D: x, 2D: y, 3D: z) | pairs |

Example — reset `disx` to the x-coordinate:

```
control_reset_dof             0  -disx
control_reset_value_linear    0  1.0 0.0 0.0
```

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
