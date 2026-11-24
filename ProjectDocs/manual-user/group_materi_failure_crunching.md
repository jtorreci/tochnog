# group_materi_failure_crunching

## Description

`group_materi_failure_crunching` (manual Professional 6.667) triggers the
element failure (slow deletion) when the **compression strain** in the
element exceeds a threshold:

> If the compression strain in an element exceeds `threshold`, the element
> is considered to be fail. The element will be slowly deleted. It is
> totally deleted if the `delete_time` has passed.

The Professional spelling includes the "n" (`crunching`). The GNU fork had
the typo `group_materi_failure_cruching` (missing "n"), which is kept as a
`db_number` alias so legacy input files keep working.

## Uso

```
group_materi_failure_crunching 0
                        0.5     ( threshold )
                        100.    ( delete_time )
```

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_failure_crunching` | `threshold delete_time` | Failure threshold on the compression strain and the deletion time (the element is slowly deleted over `delete_time`, then fully removed). |

## Física y convención de signos

The criterion is evaluated on the **most compressive principal strain**
(`tmp = min eigenvalue`, negative in compression). The element fails when
`tmp > threshold`:

- with a POSITIVE `threshold` the element never fails (`tmp <= 0`);
- with a NEGATIVE `threshold` the failure triggers as soon as the
  compression strain is above `|threshold|` — and because `tmp = 0` at
  zero strain, a negative threshold fires at the FIRST timestep even
  before any load (faithful GNU behavior, unchanged from the jan-2014
  fork; the threshold is meant to be given as a negative number of the
  "crushing strain").

The failing element is recorded in `element_delete_times` and deleted
gradually (`element_delete_factor` ramps down over `delete_time`).

## Validation

- `mcrunch.dat` / `mcrunch_low.dat` (validation-suite/test-2014): confined
  compression (oedometer), `E = 1000`, `nu = 0.3`, `eps_zz = -0.001`.
  - threshold `0.5`: no failure, elastic response intact,
    `sigma_xx = -0.5769` EXACT;
  - threshold `-0.5`: the element is marked for deletion at the first
    timestep (`element_delete_times[0] = 0.05`) — A/B discrimination.

## Notas

- Requires `materi_strain_total` in the initialization part.
- `group_materi_failure_cruching` (GNU typo) still works via the alias.
