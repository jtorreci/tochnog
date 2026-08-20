# groundflow_consolidation_apply

## Description

Global switch that enables or disables the consolidation term in the
groundflow equation. Consolidation couples the material divergence (volume
change of the solid skeleton, from `materi_velocity`) into the pore-pressure
equation.

- `-yes` (default): the material divergence part in the groundflow equation is
  included — the standard behaviour for geotechnical consolidation.
- `-no`: the material divergence part is skipped; the flow equation runs
  without the consolidation coupling.

The record is global (`no_index`). Use `control_groundflow_consolidation_apply`
to switch per timestep, or `group_groundflow_consolidation_apply` per element
group.

## Usage

```
groundflow_consolidation_apply <switch>
```

`switch` is `-yes` or `-no`. Default is `-yes`.

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `switch`  | `-yes` include consolidation, `-no` skip it. Default `-yes`.   |

## Example

```
groundflow_consolidation_apply -no
```

For a geotechnical calculation where consolidation is not wanted (e.g. purely
hydraulic flow with a prescribed skeleton velocity), the material divergence
term is omitted and the pore pressures follow the pure flow equation.
