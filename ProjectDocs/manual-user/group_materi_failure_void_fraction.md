# group_materi_failure_void_fraction

## Description

`group_materi_failure_void_fraction` (manual Professional 6.671) triggers
the element failure (slow deletion) when the **void fraction** in the
element exceeds a threshold:

> If the void fraction in an element exceeds `threshold`, the element is
> considered to be fail. The element will be slowly deleted. It is totally
> deleted if the `delete_time` has passed.

The Professional spelling uses the underscore (`void_fraction`). The GNU
fork had `group_materi_failure_voidfraction` (no underscore), kept as a
`db_number` alias for legacy input files.

## Uso

```
group_materi_failure_void_fraction 0
                        0.9     ( threshold )
                        100.    ( delete_time )
```

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_failure_void_fraction` | `threshold delete_time` | Failure threshold on the void fraction dof and the deletion time. |

## Física

The criterion compares `|void_fraction|` (the `materi_void_fraction` dof)
against `threshold`. The void fraction is an unknown of the calculation
(declared in the initialization part) that evolves with the plastic
volumetric strain; its initial value is set with `node_dof`. The failing
element is recorded in `element_delete_times` and slowly deleted over
`delete_time`.

## Validation

- `mvoid.dat` / `mvoid_low.dat` (validation-suite/test-2014): confined
  compression (oedometer) with initial void fraction `0.3`
  (`node_dof -ra -from 1 -to 4 -ra`).
  - threshold `0.9`: `0.3 > 0.9` false, element intact,
    `sigma_xx = -0.5769` EXACT;
  - threshold `0.2`: `0.3 > 0.2` true at the first timestep,
    `element_delete_times[0] = 0.05` — A/B discrimination.

## Notas

- Requires `materi_void_fraction` in the initialization part.
- `group_materi_failure_voidfraction` (GNU name) still works via the alias.
