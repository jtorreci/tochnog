# control_groundflow_nonsaturated_apply

## Description

Timestep-level switch that enables or disables the non-saturated ground water
flow data (e.g. the van Genuchten model) for the timesteps with the same
index.

- `-no`: non-saturated ground flow data (e.g. van Genuchten) is not applied for
  the selected timesteps; only saturated data is used.
- `-yes`: non-saturated data is applied.

When `control_groundflow_nonsaturated_apply` is not specified for a given
timestep, the global `groundflow_nonsaturated_apply` record is used instead
(default `-yes`).

## Usage

```
control_groundflow_nonsaturated_apply <index> <switch>
```

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `index`   | Timestep (control) index to which the switch applies.         |
| `switch`  | `-yes` or `-no`.                                              |

## Example

```
control_groundflow_nonsaturated_apply 0 -no
```

The van Genuchten model is skipped for the timesteps of control index 0 (e.g.
an initial saturated phase), while later timesteps keep the non-saturated
behaviour.
