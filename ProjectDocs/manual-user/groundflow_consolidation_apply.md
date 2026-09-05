# groundflow_consolidation_apply

## Description

Global switch that enables or disables the consolidation term in the
groundflow equation. Consolidation couples the material divergence (volume
change of the solid skeleton, from `materi_velocity`) into the pore-pressure
equation.

- `-yes`: the material divergence part in the groundflow equation is included —
  the coupling needed for geotechnical consolidation analyses.
- `-no` (default): the material divergence part is skipped; the flow equation
  runs without the consolidation coupling, and the pore pressures follow the
  pure flow equation.

The record is global (`no_index`). Use `control_groundflow_consolidation_apply`
to switch per timestep, or `group_groundflow_consolidation_apply` per element
group.

The default is `-no`, matching the Professional (manual 6.556: "Default switch
is -no"). The GNU legacy default was `-yes`; it was aligned with the
Professional in the u-p consolidation sprint (2026-09-05, commit `0560dcf`)
because it made every materi + groundflow model without an explicit switch run
a spurious consolidation transient (measured on ground14/15/16 of the corpus:
the Professional binary, without the record, reaches the drained steady state
inside the 1 s window of the test, while the legacy GNU default still carried a
large excess pore pressure at t = 1 s).

## Usage

```
groundflow_consolidation_apply <switch>
```

`switch` is `-yes` or `-no`. Default is `-no`.

## Parameters

| Parameter | Meaning                                                      |
|-----------|--------------------------------------------------------------|
| `switch`  | `-yes` include consolidation, `-no` skip it. Default `-no`.  |

## Example

```
groundflow_consolidation_apply -yes
```

For a real consolidation analysis (deforming skeleton coupled to the pore
pressure, e.g. an oedometer or a foundation settlement), set the switch
explicitly to `-yes`.

## Notes

- The switch only matters when the model combines `materi_velocity` (or the
  displacement/velocity dofs of the skeleton) with `groundflow_pressure`.
  Purely hydraulic models are unaffected.
- In a transient consolidation analysis also set `inertia_apply -yes` when the
  groundflow capacity storage term `C·ṗ` must be active (the same convention
  as the Professional; the corpus consolidation examples `large2`/`large3` set
  both records).
