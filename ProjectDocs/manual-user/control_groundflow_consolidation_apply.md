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
(default `-no`, matching the Professional manual 6.556 "Default switch is
-no").

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
groundflow_consolidation_apply -yes
control_groundflow_consolidation_apply 0 -no
```

Consolidation is active in general (the global `-yes`), but skipped for the
timesteps of control index 0 (e.g. an initial hydraulic phase).

## Notes

- The default changed from `-yes` (GNU legacy) to `-no` (Professional) in the
  u-p consolidation sprint (2026-09-05, commit `0560dcf`): without the record,
  a materi + groundflow model no longer runs a spurious consolidation
  transient. Groundflow models that want the coupled consolidation must set
  the switch to `-yes` explicitly (like `large2`/`large3` of the corpus do).
