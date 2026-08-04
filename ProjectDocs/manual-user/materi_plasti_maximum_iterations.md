# materi_plasti_maximum_iterations

## Descripción

`materi_plasti_maximum_iterations` sets a maximum number of plastic
iterations per integration point. It overrides the built-in default (10 for
viscoplastic materials, `MAX_ITER` otherwise). It is useful to control
convergence: a low limit stops slow or non-converging plasticity loops
early, a high limit allows harder problems to converge.

## Uso

Place it in the data part, grouped by material group index:

```
materi_plasti_maximum_iterations 0 5
```

## Parámetros

| Parameter | Meaning                                                          |
|-----------|------------------------------------------------------------------|
| `0`       | Material group index the limit applies to.                       |
| `N`       | Maximum number of plastic iterations per integration point. Must be at least 1. |

## Ejemplo

Run a Von Mises plastic problem with a tight iteration budget:

```
control_geometry
   cartesian
control_time
   end 1.0
   dt 0.1
control_print
   history 0
   step 10
materi_elasti_young
   0 2.e7
materi_elasti_poisson
   0 0.3
materi_density
   0 2500.
materi_plasti_vonmis
   0 1.e3 0.
materi_velocity
   0 -yes
materi_plasti_maximum_iterations
   0 5
```
