# group_materi_elasti_camclay_pressure_min

## Description

`group_materi_elasti_camclay_pressure_min index pressure_min` (manual
Professional 6.647) specifies a minimal allowed value for the pressure
in the calculation of the bulk modulus for the camclay model:

    K = (1 + e) * p / kappa

with p the current (compression-positive) mean pressure. In the
calculation, pressures below `pressure_min` are set to `pressure_min`.
This prevents numerical problems for very low bulk modulus K values
(the unclamped K -> 0 as p -> 0, and becomes negative for tensile
states). The index specifies the element_group, see `element_group`.

The record is a modifier of the camclay elastic laws
`group_materi_elasti_camclay_g` (6.645) and
`group_materi_elasti_camclay_poisson` (6.646).

## Uso

In the data part:

```
group_materi_elasti_camclay_g 0  1000.
group_materi_elasti_camclay_pressure_min 0  10.
group_materi_plasti_camclay 0
                        0.8    ( m )
                        0.02   ( kappa )
                        0.2    ( lambda )
                        3.     ( N )
```

## Parámetros

| # | Meaning |
|---|---------|
| pressure_min | minimal allowed pressure for the camclay bulk modulus (pressures below are clamped up to it) |

## Validation

- `mc_pressure_min` (confined compression with camclay elastic, initial
  p = 0.001, pressure_min = 10): K clamped to (1+e)*10/kappa = 750;
  one step gives `sigxx = -0.0093333` and `sigyy = -0.2093333` EXACT
  (analytic: dp = K*tr(eps) = 0.075, sigxx = -0.001 + lambda*tr,
  sigyy = -0.001 + lambda*tr + 2G*eps_yy).
- A/B `mc_pressure_min_off` (no record): K = (1+e)*0.001/kappa = 0.075
  (near-singular), the response degenerates: `sigxx = +0.06566`
  (tension, sign flip) — the numerical problem the manual prevents.
