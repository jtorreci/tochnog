# check_target

## Descripción

`check_target` controls how the `target_item` / `target_value` records are
handled at the end of the calculation. By default (`-yes`) an unmet target
aborts the run with an error, as usual. With `check_target -no` unmet targets
no longer abort: they are downgraded to a "Note" entry in `tn.log` and the run
continues. This is useful while probing variations of the input file, when you
want the analysis to finish even though the target is not met.

## Uso

Place it in the data part, as a keyword line:

```
check_target -no
```

## Parámetros

| Parameter | Meaning                                                                 |
|-----------|-------------------------------------------------------------------------|
| `-yes`    | Unmet targets abort the run with an error (default).                   |
| `-no`     | Unmet targets are logged as a Note in `tn.log` and the run continues. |

## Ejemplo

Minimal input that keeps running even though its target is not met:

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
check_target
   -no
```

When a `target_item` / `target_value` record is not satisfied, `tn.log`
contains a line similar to:

```
Note in calculation with data file <file>.
Target value (check_target -no) not met for: <data_item> index <index> (wanted <value>).
```

and the run completes instead of aborting.
