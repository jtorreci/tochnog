# check_used

## Descripción

`check_used` reports the data items that were defined in the data file but
were never read during the calculation. This helps detect orphaned data:
typos in keyword names, conditions that never apply, or leftover entries
from an older input file. Enabling it makes the run fail-fast on mistakes
that would otherwise pass silently.

## Uso

Place it in the data part, as a keyword line:

```
check_used -yes
```

## Parámetros

| Parameter | Meaning                                                    |
|-----------|------------------------------------------------------------|
| `-yes`    | Enable the report. At the end of the run a list of unused data items is printed to the standard output. |
| `-no`     | Disable the report (default).                              |

## Ejemplo

Minimal input that defines an unused item on purpose:

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
materi_velocity
   0 -yes
# unused on purpose: material is never used by any element
materi_plasti_mohrcoul
   0 1.e3 0.3 0.2
check_used
   -yes
```

Run it and the log will list `materi_plasti_mohrcoul` as defined but not used.
