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

- In Professional the limit defaults to 0 even when the record is absent.
  In the GNU the record is REQUIRED to activate the limit: without it
  nothing is clamped. This keeps backwards compatibility with existing
  GNU models that legitimately solve positive pressures (e.g. the
  seepage example from the Professional manual itself, which develops
  positive heads).
- Nodes with a prescribed pressure (`bounda_dof`/`bounda_unknown`, or the
  phreatic `_static` condition) keep their prescribed value; only solved
  (free) nodes are cut. The manual text "any higher value resulting from
  the groundflow equations" is read literally.

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
