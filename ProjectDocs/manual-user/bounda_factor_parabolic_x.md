# bounda_factor_parabolic_x

## Description

Multiplies the load of a boundary by a factor that is quadratic in the `x`
coordinate of the node:

```
factor = a0 + a1*x + a2*x^2
```

The load applied to each node is multiplied by this factor, allowing parabolic
spatial load distributions along `x` (e.g. a parabolic pressure profile).
Requires `bounda_time` for the same boundary (`iboun`).

## Usage

```
bounda_factor_parabolic_x <iboun> <a0> <a1> <a2>
```

## Parameters

| Parameter | Meaning                                      |
|-----------|----------------------------------------------|
| `iboun`   | Boundary number to which the factor applies. |
| `a0`      | Constant term of the factor.                 |
| `a1`      | Coefficient of `x`.                          |
| `a2`      | Coefficient of `x^2`.                        |

## Example

Apply a parabolic factor `1. + x + x^2` to the load of boundary 1:

```
bounda_time                 1  0. 10. 100. 10.
bounda_factor_parabolic_x   1  1. 1. 1.
```

At a node with `x=1.` the factor is `3.` and the load applied is `3*10. = 30.`;
at `x=2.` the factor is `7.`.
