# bounda_time_until_force

## Descripción

Limits the reaction force of a node with a prescribed velocity (or
displacement). If the reaction force from the previous step exceeds
`max_force`, the prescribed velocity of the current step is reduced so that
the applied force stays controlled without exceeding the limit.

When active, the prescribed velocity is multiplied by:

```
factor * ( max_force / reaction )
```

Useful to apply a force-controlled load without exceeding a limit while
keeping the `bounda_unknown`/velocity formulation. Requires `bounda_time`
for the same boundary (`iboun`), which provides the target load amplitude.

## Uso

```
bounda_time_until_force <iboun> <max_force> <factor>
```

## Parámetros

| Parámetro   | Significado                                                         |
|-------------|---------------------------------------------------------------------|
| `iboun`     | Boundary number of the prescribed-velocity record to limit.         |
| `max_force` | Maximum allowed reaction force.                                     |
| `factor`    | Scaling factor applied to the reduction, must be `0 < factor <= 1`. |

Only affects prescribed velocity/displacement boundaries (`bounda_unknown`),
not applied forces (`bounda_force`).

## Ejemplo

Limit the reaction of boundary 1 to `5.0`, reducing the prescribed velocity
by half of the required factor:

```
bounda_unknown  1  -geometry_point 1 -velx
bounda_time     1  0. 0.1 100. 0.1
bounda_time_until_force  1  5.0 0.5
```

The reaction force of boundary 1 never exceeds `5.0`; when the raw reaction
would exceed it, the velocity is scaled down by `0.5 * (5.0 / reaction)`.
