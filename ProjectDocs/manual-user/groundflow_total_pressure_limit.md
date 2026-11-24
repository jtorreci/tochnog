# groundflow_total_pressure_limit

## Description

Maximum allowed total pressure (manual Professional 6.587): any higher
calculated total pressure is cut off to this value. **Default the limit is
0**, so the record is only needed to change it.

What the limit does:

- The **calculated total pressure** (`post_calcul -groundflow_pressure
  -total_pressure`, `-to_pres`) is capped at the limit at every node,
  prescribed nodes included. The same cap applies wherever the program
  uses the total pressure (pore pressure on the skeleton, wet/dry density,
  `-total_pressure` boundary conversions).
- The **groundflow velocity** is 0 where the total pressure at an
  integration point (interpolated from the capped nodal total pressures)
  reaches the limit: there is no water there (the dry rule of manual
  2.4.3). With the default 0 this means a zero Darcy velocity wherever the
  total pressure is 0 or positive (at and above the water table).
- The **`pres` dof is never cut**: the groundflow equation, including the
  consolidation term, is solved unchanged. The limit only acts on what is
  derived from the total pressure.

Typical use: the default 0 models soil that cannot carry suction (positive
total pressure in the tension-positive convention). Set a large value (for
example `1.e10`) to keep positive total pressures and the Darcy velocity
everywhere; set a negative value to cap the total pressure below 0.

Sign convention: tension positive. With `force_gravity (0,-1)` the pore
pressure below the water level is negative (compression), so the limit
caps the positive (suction) side; compression passes through.

`groundflow_pressure` must be an unknown of the analysis.

## Usage

```
groundflow_total_pressure_limit <limit>
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `limit`   | Maximum allowed calculated total pressure. Default 0 (record absent). |

## Example

One quad4 with a prescribed solid velocity (divergence 2), consolidation
and capacity storage active, pres 0 prescribed at the bottom; the free top
nodes develop `pres` 1.6667 (validation tests
`test-2026/groundflow_total_pressure_limit_default1`, `_high1`,
`_negative1`, `test-2014/groundflow_total_pressure_limit`, `_dry`):

```
groundflow_consolidation_apply   -yes
inertia_apply   -yes
groundflow_total_pressure_limit  0.5
post_calcul -groundflow_pressure -total_pressure
```

| input | free-node `pres` | total pressure top / bottom | `gvely` |
|---|---|---|---|
| no record (limit 0) | 1.6667 | 0 / 0 | 0 |
| limit 0 | 1.6667 | 0 / 0 | 0 |
| limit 0.5 | 1.6667 | 0.5 / 0 | 0.16667 |
| limit 100 | 1.6667 | 1.6667 / 0 | 0.16667 |
| limit -0.5 | 1.6667 | -0.5 / -0.5 | 0 |

Corpus `ground1` (pres 1 at the bottom, 0 at the top, no gravity): the
total pressure is 1 / 0.5 / 0, capped to 0 at every node, so the limit is
reached everywhere and `gvely` is 0; with `groundflow_total_pressure_limit
1.e10` it is -0.05 (`test-2026/groundflow_velocity_dry1`, `_dry2`).

## Differences with the Professional version

None known. The behaviour above was measured on the Professional
02-08-2026 (F4 T1 probes: the five rows of the table, positive and
negative pressure fields, and `ground1`) and the GNU gives the same
values.

History: until F4 (2026-10-09) the GNU needed the record to activate the
limit, cut the solved `pres` dof of the free nodes, skipped the
consolidation term when the limit was 0 and the pressure 0, and capped the
post-processed total pressure with `groundflow_pressure_atmospheric`
instead. Aligned with manual 6.587 and the Professional by project-owner
decision 2026-10-07 (former DIVERGENCES TPL and K10). Note that manual
2.4.3 says the consolidation part is skipped when the limit is 0 and the
total pressure is 0; the Professional does not skip it (free-node `pres`
1.6667 with the default limit), and the GNU follows the Professional.

- **Deliberate divergence from the manual (decision 2026-10-09, the "ley seca"):** the dry rule of manual 2.4.3 (limit 0 and pres 0 skip the consolidation term) is not applied, as in the Professional 02-08-2026. See `ProjectDocs/DIVERGENCES-PROFESSIONAL.md` (DRY-RULE).
