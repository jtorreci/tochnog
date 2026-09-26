# control_plasti_apply

## Description

`control_plasti_apply index switch` (manual Professional 6.257): with
`-no`, any plasticity data is ignored for the timesteps of the
`control_timestep` record with the same index. It is the per-timestep form
of `plasti_apply` (6.899); the scope of "plasticity data" is described in
`plasti_apply.md`.

## Usage

```
control_plasti_apply index switch
```

| Parameter | Meaning |
|-----------|---------|
| `index` | the `control_timestep` index |
| `switch` | `-yes` (default) / `-no` |

The record overrules a global `plasti_apply` for its index; other indices
use `plasti_apply` (default `-yes`). Plastic strains of earlier timesteps
are kept.

## Example

```
group_materi_plasti_tension_direct 0 1.
control_timestep 0 0.1 0.5
control_timestep 1 0.1 0.5
control_plasti_apply 1 -no
```

Tension cut-off at 1 in timestep 0, elastic increment in timestep 1:
sigyy 11.989 (Professional 02-08-2026 and GNU,
`validation-suite/test-2026/control_plasti_apply1.dat`).
`control_plasti_apply2.dat` gives the same result with `plasti_apply -no`
overruled by `control_plasti_apply 0 -yes`.

## Differences with the Professional version

See `plasti_apply.md`.
