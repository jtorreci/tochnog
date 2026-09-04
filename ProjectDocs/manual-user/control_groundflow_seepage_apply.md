# control_groundflow_seepage_apply

## Description

`control_groundflow_seepage_apply` (manual Professional 6.105) is the
per-control-step switch of the groundflow seepage faces:

```
control_groundflow_seepage_apply index switch
```

Set the switch to `-no` to disable the seepage faces during the
selected control steps (e.g. a linear elastic stage before the
nonlinear seepage stage), `-yes` to enable them.

## Uso

```
control_timestep                      10  1. 1.
control_groundflow_nonsaturated_apply 10  -no
control_groundflow_seepage_apply      10  -no
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `index` | Control step index (must match a `control_timestep`). |
| `switch` | `-yes` (default): seepage faces active; `-no`: disabled. |

## Verification

Registered and validated for parsing (tutorial_2, which uses it, gets
past the record). Requires a `groundflow_pressure` analysis (checked).
NOTE: the gate is not yet consumed by the seepage machinery in
`bounda.cc`; the tutorial additionally needs its `include`d mesh which
the corpus harness does not stage.
