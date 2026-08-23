# groundflow_pressure_atmospheric

## Description

Upper cap for the static and total pore pressure derived from the phreatic
level in `groundflow_phreatic_coord()`. Any value more positive than the
cap is truncated to the cap. It is a cap, not a gauge conversion: nothing
is subtracted from the pressure, values are simply clipped.

Sign convention: with gravity pointing down (`force_gravity 0. -1.`), the
static pore pressure below the phreatic level is negative (compression) and
suction above the level is positive. With the default cap 0, suction above
the phreatic surface is removed (the classic "no suction" behaviour) while
compression below the level passes through untouched. A positive cap keeps
suction up to that value.

This keyword is inherited from the GNU codebase (it was already present in
the sfnet 2014 sources); it does not exist in Tochnog Professional 2024.

Note: `bounda_water` computes the hydrostatic pore pressure directly and
deliberately bypasses this cap (see `manual-developer/bounda_water.md`).

## Usage

```
groundflow_pressure_atmospheric <pa>
```

## Parameters

| Parameter | Meaning                                                       |
|-----------|---------------------------------------------------------------|
| `pa`      | Cap on the static/total pore pressure. Default `0`.           |

## Example

```
force_gravity           0. -1.
groundflow_density      1.
groundflow_pressure_atmospheric  0.5
groundflow_phreatic_level_multiple 0  1.
groundflow_phreatic_level_multiple_element_group 0  0
groundflow_phreatic_level_multiple_static 0  -yes
```

Column with nodes at y=0,1,2 and phreatic level y=1: with `_static -yes`
the imposed pressure at y=0 is the full compression `-1.0` (negative, not
affected by the cap) and at y=2 the suction `+1.0` is capped to `+0.5`.
Without the keyword (default 0) the suction node gets `0.0`.
