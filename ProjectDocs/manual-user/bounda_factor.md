# bounda_factor

## Description

Multiplies the load of a boundary by a factor that depends on the spatial
coordinates of the node. The factor is linear in the coordinates:

```
1D:  a0 + a1*x
2D:  a0 + a1*x + a2*y
3D:  a0 + a1*x + a2*y + a3*z
```

Only the coefficients relevant to the current number of dimensions are used.
The load applied to each node is multiplied by this factor, allowing spatially
varying boundary loads (e.g. a triangular pressure distribution). Requires
`bounda_time` for the same boundary (`iboun`).

## Usage

```
bounda_factor <iboun> <a0> <a1> [<a2> [<a3>]]
```

## Parameters

| Parameter | Meaning                                     |
|-----------|---------------------------------------------|
| `iboun`   | Boundary number to which the factor applies.|
| `a0`      | Constant term of the factor.                |
| `a1`      | Coefficient of `x`.                         |
| `a2`      | Coefficient of `y` (2D/3D).                 |
| `a3`      | Coefficient of `z` (3D).                    |

## Example

Apply a linear ramp in `x` to the load of boundary 1 (factor `1. + 2.*x`):

```
bounda_time   1  0. 10. 100. 10.
bounda_factor 1  1. 2.
```

At a node with `x=1.` the factor is `3.` and the load applied is `3*10. = 30.`;
at `x=2.` the factor is `5.`.
