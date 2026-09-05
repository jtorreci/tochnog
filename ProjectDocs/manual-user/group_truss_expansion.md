# group_truss_expansion

## Description

With `group_truss_expansion` you can specify the thermal expansion
coefficient of the truss material (manual Professional 6.775). A
temperature increment `dT` produces a free thermal incremental length of
`alpha * dT * initial_length`; the truss force builds up only from the
mechanical stretch that remains after subtracting the free thermal
expansion. A truss that is prevented from expanding while its
temperature rises therefore develops a compressive force of
`-E*A*alpha*dT`.

The temperature of the truss is the average of its nodal temperatures
(the `condif_temperature` dof when a temperature analysis is active).

## Usage

```
group_type <group_index> -truss
group_truss_expansion <group_index> <alpha>
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `alpha`   | Thermal expansion coefficient (1/K); a temperature increment `dT` gives a free thermal incremental length `alpha*dT*L0`. |

## Example

See `truss11.dat`: one truss of E=A=1 and alpha=1 heated to T=2 while
its ends are fixed develops `element_truss_force = -2` (verified against
the Professional `.dbs`).
