# groundflow_nonsaturated_apply

## Description

Global switch that enables or disables the non-saturated ground water flow
data (e.g. the van Genuchten model). When `-no`, only the saturated data
(`group_groundflow_capacity`, `group_groundflow_permeability`) is used.

- `-yes` (default): non-saturated data (e.g. van Genuchten) is applied.
- `-no`: only saturated data is used.

This is done for all timesteps. Use `control_groundflow_nonsaturated_apply`
to switch per timestep.

## Usage

```
groundflow_nonsaturated_apply <switch>
```

`switch` is `-yes` or `-no`. Default is `-yes`.

## Parameters

| Parameter | Meaning                                                          |
|-----------|------------------------------------------------------------------|
| `switch`  | `-yes` apply non-saturated data, `-no` saturated only. Default `-yes`. |

## Example

```
groundflow_nonsaturated_apply -no
```

The van Genuchten model is ignored and the calculation runs saturated.
