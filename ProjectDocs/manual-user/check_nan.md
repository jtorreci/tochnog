# check_nan

## Descripción

`check_nan` inspects the `node_dof` results for NaN (Not A Number) values. If a
NaN is detected, the run aborts indicating the node and the dof. A NaN in the
solution usually means the analysis diverged or there is an error in the input.
It is a diagnostic tool, disabled by default.

## Uso

Place it in the data part, as a keyword line:

```
check_nan -yes
```

## Parámetros

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `-yes`    | Enable the NaN check. Aborts on the first NaN found in `node_dof`. |
| `-no`     | Disable the check (default).                                   |

## Ejemplo

Minimal input that enables the NaN check:

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
check_nan
   -yes
```

The check runs at the end of every step close. If a NaN is found the run stops
with an error like:

```
Error: NAN detected in node_dof of node <inod> dof <iuknwn>.
The solution may have diverged.
```
