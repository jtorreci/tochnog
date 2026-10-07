# groundflow_total_pressure_limit

## Description

Maximum allowed pore (total) pressure value. Any higher value resulting
from the groundflow equations is cut off to this value. In addition, when
the limit is 0 and the pressure of a node is 0, Tochnog assumes there is
no water there and the consolidation term (material divergence) is skipped
for the dry elements (manual Professional, section 2.4.3): water cannot
sustain suction, so a dry zone does not couple to the flow equations.

Typical use: dewatering/excavation models where the water table drops and
parts of the mesh become dry. With `limit 0` the dry nodes stay at 0
instead of developing unphysical suction.

Sign convention: with `force_gravity (0,-1)` compression below the water
level is negative and suction is positive, so the limit caps the
(positive) suction side; negative compression passes through.

`groundflow_pressure` must be an unknown of the analysis.

### Differences with the Professional version

Measured 2026-10-07 on the Professional 02-08-2026 (probe: one quad4,
prescribed solid velocity with divergence 2, `groundflow_consolidation_apply
-yes`, `inertia_apply -yes`, pres 0 prescribed at the bottom, free top
nodes; uncapped free-node pres 1.6667). The Professional manual (6.587)
says: "the maximum allowed total pressure value. Any higher value resulting
from the groundflow equations will be cutoff to this value. Default the
limit is set to 0."

| input | GNU free-node `pres` | Pro free-node `pres` | GNU total pressure (top / bottom) | Pro total pressure (top / bottom) | Pro `gvely` |
|---|---|---|---|---|---|
| no record | 1.6667 | 1.6667 | 0 / 0 | 0 / 0 | 0 |
| limit 0 | 0 | 1.6667 | 0 / 0 | 0 / 0 | 0 |
| limit 0.5 | 0.5 | 1.6667 | 0 / 0 | 0.5 / 0 | 0.16667 |
| limit 100 | 1.6667 | 1.6667 | 0 / 0 | 1.6667 / 0 | 0.16667 |
| limit -0.5 | -0.5 | 1.6667 | -0.5 / 0 | -0.5 / -0.5 | 0 |

- **What is cut and the default** (known difference, pending project-owner
  decision; DIVERGENCES TPL). The Professional reads the record (it is not
  a GNU-only record). It caps the CALCULATED total pressure
  (`post_calcul -groundflow_pressure -total_pressure`, which the
  Professional requires whenever a pressure dof exists) at every node,
  prescribed nodes included, and never cuts the `pres` dof itself; the
  default limit is 0, so a positive total pressure is reported as 0 even
  without the record. Where the total pressure reaches the limit the
  Professional reports a zero groundflow velocity (the dry rule of manual
  2.4.3); this is the cause of the `ground1` velocity difference (AUDIT
  K10): with `groundflow_total_pressure_limit 1.e10` the Professional
  gives -0.05 like the GNU. With limit 0 the Professional still develops
  the free-node `pres` (1.6667).
  The GNU needs the record to activate the limit (no record: nothing is
  cut, kept for existing GNU models with positive heads), cuts the
  SOLVED `pres` dof of the free nodes only, and with limit 0 + pres 0
  skips the consolidation term (pres stays 0). Its post-processed total
  pressure does not depend on this record: it is capped by
  `groundflow_pressure_atmospheric` (default 0). For negative pressures
  (below the limit) both programs agree (-1.6667).
  Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).

## Usage

```
groundflow_total_pressure_limit <limit>
```

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `limit`   | Maximum allowed solved pore pressure. Use 0 for the dry-soil behaviour. |

## Example

```
groundflow_total_pressure_limit  0.
```

With a velocity-divergence source that would drive the free-node pressure
to ~1.67, the free nodes stay at 0 (dry, consolidation term skipped) and
nodes with prescribed pressure keep their value.
