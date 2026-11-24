# post_point_move

## Descripción

`post_point_move` makes a post point follow the material particle: instead
of staying at its initial coordinate, its position is updated every time
step with the interpolated velocity times the time step. It is useful to
track the history of a material point in large-deformation problems.

## Uso

Place it in the data part:

```
post_point_move -yes
```

Requires `materi_velocity` to be active in the model.

## Parámetros

| Parameter | Meaning                                                          |
|-----------|------------------------------------------------------------------|
| `-yes`    | Update each post point coordinate each step with the interpolated velocity field. |
| `-no`     | Keep post points fixed at their initial coordinates (default).  |

## Ejemplo

Track the velocity/stress history of a material particle in a large-slip run:

```
control_geometry
   cartesian
control_time
   end 10.0
   dt 0.1
control_print
   history 0
   step 10
post_point
   0 0.5 0.5
post_point_move
   -yes
materi_elasti_young
   0 2.e7
materi_elasti_poisson
   0 0.3
materi_density
   0 2500.
materi_velocity
   0 -yes
```
