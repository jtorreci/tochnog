# volume_factor_x

## Description

An in x-direction CHANGING volume factor for elements (manual
Professional 6.1094): a piecewise-constant thickness (2D) or
cross-sectional area (1D) keyed on the global x coordinate. Left from
the first x the factor is 1, between x_i and x_i+1 it is the listed
factor, and right from the last x it is 1 again. It multiplies the
polynomial `volume_factor` and the per-group `group_volume_factor`
when those are also present.

## Syntax

```
volume_factor_x x0 fac01 x1 fac12 ... xn
```

An odd number of values: positions and factors alternating. Example:
`volume_factor_x 4. 2. 8.` gives factor 1 up to x=4, factor 2 from
x=4 to x=8, factor 1 beyond x=8.

## Example

The `tslv_vfx` test of the validation suite: an 8x1 quad4 chain
(E=1000, total pull N=2e-2) with `volume_factor_x 4. 2. 8.` - the
right half carries thickness 2. Series statics: sigxx = 2e-2 in the
left half, 1e-2 in the right half, tip displacement
N*4/(E*1) + N*4/(E*2) = 1.2e-4 (measured 1.15e-3 velocity at dt=0.1,
within the Poisson band; without the record 1.59e-3, the A/B).
