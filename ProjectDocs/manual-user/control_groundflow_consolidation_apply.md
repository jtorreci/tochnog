# control_groundflow_consolidation_apply

## Description

Timestep-level switch that enables or disables the consolidation term (the
material divergence part of the groundflow equation) for the timesteps with
the same index.

- `-no`: the material divergence part in the groundflow equation is skipped for
  the selected timesteps.
- `-yes`: consolidation is active.

When `control_groundflow_consolidation_apply` is not specified for a given
timestep, the global `groundflow_consolidation_apply` record is used instead
(default `-yes`).

## Usage

```
control_groundflow_consolidation_apply <index> <switch>
```

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `index`   | Timestep (control) index to which the switch applies.         |
| `switch`  | `-yes` or `-no`.                                              |

## Example

```
control_groundflow_consolidation_apply 0 -no
```

Consolidation is skipped for the timesteps of control index 0 (e.g. an initial
hydraulic phase), while later timesteps keep the default behaviour.
