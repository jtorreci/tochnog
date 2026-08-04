# bounda_time_on_off

## Descripción

Applies the values of `bounda_time` periodically (cyclic loads). The load is
active only during the "on" window of every cycle:

```
fmod( time_total, time_on + time_off ) < time_on
```

During the "off" window the load is 0. Requires `bounda_time` for the same
boundary (`iboun`), which provides the load amplitude; this keyword only
switches that amplitude on and off.

## Uso

```
bounda_time_on_off <iboun> <time_on> <time_off>
```

## Parámetros

| Parámetro  | Significado                                                        |
|------------|--------------------------------------------------------------------|
| `iboun`    | Boundary number to which the `bounda_time` record applies.         |
| `time_on`  | Duration of the active (loaded) window, must be > 0.               |
| `time_off` | Duration of the inactive (zero load) window, must be >= 0.         |

The cycle period is `time_on + time_off`.

## Ejemplo

Apply the load defined by `bounda_time 1` only 40% of each 1.0 time-unit cycle:

```
bounda_time     1  0. 100. 100. 100.
bounda_time_on_off  1  0.4 0.6
```

The load of boundary 1 is active for `0.4` time units and zero for `0.6`,
repeating periodically. Useful for cyclic loading (waves, traffic, etc.).
