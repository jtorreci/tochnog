# Integration quadrature and nodal stress recovery (Professional alignment)

Reference for **which quadrature each operation uses** in Tochnog
Professional and in this GNU version, and for **how the Professional
recovers the nodal stresses** from the integration-point stresses. This is
the developer reference behind the 2026-09-13 alignment (branch
`quad/pro-integration-recovery`).

The task that produced this document: map the quadrature PER OPERATION
(not a single global default) and determine whether GNU can match the
Professional by choosing the method per operation.

## 1. Per-operation quadrature map

The key finding: the quadrature is a property of the **element integral**
(stiffness, mass, internal force, body load), and the **nodal stress
recovery is a separate operator** that is NOT tied to the quadrature. The
Professional uses interior Gauss for the element integral and a
least-squares extrapolation for the recovery; GNU historically conflated
the two by using the shape function `h` of the same quadrature for the
recovery.

| operation | Professional (binary 2-8-2026) | GNU (before) | GNU (aligned) | evidence |
|---|---|---|---|---|
| element stiffness, quad4 / hex8 | 2x2 / 2x2x2 **Gauss** | 2x2 / 2x2x2 Lobatto (`-minimal` default) | Gauss | static probes, 12 digits (prior sprint) |
| element stiffness, quad9 / hex27 | 3x3 / 3x3x3 **Gauss** | Lobatto | Gauss | `.dbs element_intpnt_method -gauss`, `npoint` 9/27 |
| mass matrix | **row-sum lumped** (quadrature-invariant) | row-sum lumped | unchanged | distorted-element dynamic probe |
| internal force vector | same rule as the stiffness (one `pol()` rule) | Lobatto | Gauss | — |
| body / gravity load | same rule | Lobatto | Gauss | — |
| edge / surface load | boundary integration; interface elements keep nodal (Lobatto) | boundary integration | unchanged | GNU `area.cc`; interfaces do not call `pol()` |
| stress recovery (IP -> node) | **per-element least-squares fit + arithmetic mean** (section 2) | `h`-weighted lumped average | Professional operator | `.dbs element_intpnt_dof -> element_dof -> node_dof` |
| post-processing / nodal averaging | arithmetic mean over the adjacent elements (`node_lhside` = element count) | volume-weighted | arithmetic mean | `.dbs node_lhside` |

Professional manual 6.623 says `-minimal`/`-maximal` use Gauss and "in case
inertia is applied then lobatto". The 2-8-2026 binary contradicts this for
the STIFFNESS with inertia (dynamic8 matches Gauss to 6e-11). The clause
most plausibly describes the MASS, and the row-sum lumped mass is
quadrature-invariant, so GNU needs no change there.

## 2. The Professional's nodal stress recovery

The Professional `.dbs` dump exposes everything needed to identify the
operator:

- `element_intpnt_dof` — the stress at each integration point;
- `element_dof` — the per-element nodal value;
- `node_dof` — the assembled nodal value;
- `element_intpnt_h` — the shape functions at the integration points;
- `node_lhside` — the assembled weight.

Measured on the Professional's own `beam2d_1.dbs` (quad9, 45 elements),
`q4b.dbs` (quad4) and `static_hex8cube.dbs` (hex8):

| element | points | operator | residual vs `element_dof` |
|---|---|---|---|
| quad4 | 4 | **bilinear interpolation** `1,x,y,xy` (4 terms) | 1.4e-7 |
| quad9 | 9 | **uniform-weight linear least squares** `1,x,y` (3 terms) | 1.2e-12 |
| hex8 | 8 | **uniform-weight linear least squares** `1,x,y,z` (4 terms) | 1.0e-8 |
| hex27 | 27 | uniform-weight linear least squares `1,x,y,z` (assumed) | not independently verified |

The operator is FIXED in the element natural coordinates (the same 9x9
matrix fits all 45 beam2d_1 elements of different geometry to 1e-9). The
least-squares fit uses UNIFORM weights (the Gauss weights give a different
operator). The order rule is not a single clean formula: quad4 uses the
full 4-term interpolation while quad9 and hex8 use the linear fit. The
operators above were determined empirically and are the discriminating
data; hex27/quad16 were not independently exercised.

Assembly: `node_dof` = **arithmetic mean** of the adjacent elements'
`element_dof` (e.g. beam2d_1 node 3: `node_lhside` = 2 and
`node_dof = (54.022562 + 52.586786)/2 = 53.304674`). The Professional's
per-element value has `sum_points F = 1`, so each element contributes
weight 1.

## 3. GNU implementation

- `miscel.cc` — `sri_stress_recovery_weight(sri, name, nnol, inol, npoint,
  ipoint, h_inol, volume)`. Two paths:
  - `sri == 1`: the GNU-specific SRI feature, unchanged (Lagrange
    extrapolation of the quad4/hex8 normal components, volume-weighted
    assembly).
  - `sri == 0`: the Professional alignment. Returns the operator factor
    divided by the point volume. The lumped sigma-dof equation multiplies
    the weight by the volume again (`volume * weight`), so
    `volume * F/volume = F`, each element contributes `sum_points F = 1`
    and the assembly becomes the arithmetic mean. Without the volume
    cancellation the quadrature weights bias the average (measured on a
    controlled quad9 patch: the exact `sigxx = 1000` came out 1333 = 4/3).
- `polynom.cc` — `stress_recovery_gauss` (global) is 1 when the resolved
  quadrature is interior Gauss and 0 for Lobatto. `pol()` sets it; the
  recovery returns `h` for every case that is not the aligned Gauss rule.
- `general.cc` / `materi.cc` — apply the recovery to ALL stress components
  (normal and shear), exactly as the Professional's `element_dof`; the
  SRI path keeps its normal-only recovery. The `materi.cc` stress block
  was moved out of the `|h| > EPS_H` guard because the linear
  least-squares operator is dense (some `(node, point)` pairs have `h = 0`
  but a non-zero factor, e.g. the mid-edge nodes of quad9 at the centre
  Gauss point).
- `tochnog.h` / `tochnog-mod.h` — the new signature and the
  `stress_recovery_gauss` global.

### 3.1 Absolute assignment (2026-09-15)

The 2026-09-13 alignment fixed the WEIGHTS but the recovery still
ACCUMULATED the increments (`node_dof_new += node_rhside/node_lhside`), so
at a node whose elements' increments cancel the reported stress kept its
step-start value. The 2026-09-15 work unit (`nodal_stress_recovery.md`)
makes the stress dofs ASSIGN the ratio:

```c
node_dof_new[stress] = node_rhside / node_lhside;   // dof.cc
node_rhside[stress] += volume*weight*new_sig;       // materi.cc
node_lhside[stress] += volume*weight;               // materi.cc
```

A currently deleted element contributes neither term (`elem.cc` zeroes its
stress slots), matching the Professional's `node_lhside = 0` on the nodes
exclusive to a deleted element. The absolute assignment is applied only when
the stress dofs are a pure recovered field (gate in `general.cc`:
`!materi_strain_total` and the stress dof not convected).

### 3.2 Axisymmetric rule (K20, F4 T4, 2026-10-09)

`pol()` (`polynom.cc`) used to force `-GAUSS -MINIMAL` (1 point) for
every axisymmetric `materi_velocity` tensor-product element and to stop
with "Error, not available for axisymmetric analysis" for `npol > 2`
(quad9, quad16). `msf_element_rule()` (`calcul_force.cc`) replicated the
forced rule. Both blocks were removed: the axisymmetric quad4/quad9 now
use the default 2x2 / 3x3 Gauss rule (and any
`group_integration_*` override), as the Professional 02-08-2026.

Evidence (Professional 02-08-2026 probes, investigation of F4 T4):

| case | old GNU | new GNU | Professional |
|---|---|---|---|
| Akin ring `axisym2` (1 quad4): inner disx | 1.08238e-2 | 1.081776788020e-2 | 1.081776788020e-2 |
| same, inner sigzz (nodal) | 10000 (element-constant) | 10760.09 | 10759.22 |
| 4x4 quad4 ring patch, point load (`hg4`): disy | ~1e12 (hourglass, no convergence) | -4.38051133e-4 | -4.380511327e-4 |
| quad9 axisym (`axisym7`): sigxx / disx | refused | 1000 / 500.00004 | 1000 / 499.9988 |

The 1-point rule left the quad4 hourglass modes free (the plane quad4 is
saved by the 2x2 rule of the stiffness) and made the hoop strain u/r
constant over the element. The remaining `axisym2` sigzz difference is
the quad4 nodal recovery (natural-coordinate bilinear interpolation in
the GNU, physical-coordinate linear least squares in the Professional;
BACKLOG C42). `msf_axisym` (`post_calcul -materi_stress -force`): the
section moment changes from 0 to 1/960, the analytic moment of the
2*PI*r-weighted face (the build_safe tolerance 0.005 holds).
`foundation_consolidation_settlement` (Pro `other`, axisymmetric) now
runs past the point where the 1-point rule failed and needs more than
4 GB (F4 T7 memory hogs).

## 4. force5 (phreatic level)

`area.cc` looked up the phreatic level for the FIRST node of a
non-selected side before the enclosing loop's `ok` guard skipped the
side. With the Gauss solution the boundary node at x=2 deforms to
x = 2.000000000358, just outside the phreatic table [0,2], and
`table_xy` (`math.cc`) has no upper tolerance, so the run aborted
("no valid phreatic level for FORCE_ELEMENT_EDGE_WATER"). Fix: skip the
phreatic lookup for nodes that are not on the selected edge geometry
(`ok_tmp`); a non-selected side carries no water force. With the fix
force5 completes and gives `sigxx = -9.80718`, the Professional's value.

## 5. Validation (branch `quad/pro-integration-recovery`, abf9c29 + this diff)

| test | result | Professional |
|---|---|---|
| dynamic8 final disx20 | -9.98525437506e-4 | -9.985253775415e-4 (abs 6.0e-11) |
| post4 moment | 0.666667 | n/a |
| beam2d / beam2d_1 corner sigxx | 59.5766 | 59.5710 |
| force5 post_point sigxx | -9.8071801515 | -9.80718 |
| examp16 | PASS | — |
| sfnet 197 corpus | **155 PASS / 42 RUNFAIL** | baseline Lobatto 160; prior Gauss+fixA 152 |

Corpus regressions vs the Lobatto baseline (5, all outside the aligned
operations):

| test | GNU aligned | classification |
|---|---|---|
| condif10 | 17.5904 | STALE: Pro ships 18.4151 on the equivalent h-refined mesh; GNU-Gauss with h-refinement = 18.4151 |
| examp5 | 1.03522 | STALE: inverse target calibrated to the Lobatto forward value |
| ho_mech1 | 1.32827 | STALE: Pro ships the same test as `pressed_rubber_cylinder` with target 1.32791; GNU-Gauss = 1.327915 |
| strconv1 | 0.983386 | stale transient: converges to 1.0 exactly at 200 steps (20-step target calibrated to the Lobatto convection operator); no Professional reference (keyword set differs) |
| examp8 | 0.132011 | unknown: ~10% transient conduction shift, no Professional reference in the distribution |

## 6. Pending refactorings / limitations

- The recovery operators are empirically matched to the Professional
  binary; hex27 and quad16 were not independently verified (assumed linear
  least squares). A second element type would need a discriminating
  `.dbs` dump.
- The fit basis differs between quad4 (bilinear) and quad9/hex8 (linear);
  the underlying rule in the Professional was not identified, only the
  operators.
- The 3 stale targets (condif10, examp5, ho_mech1) should be recalibrated
  in the same work unit that promotes the alignment to the default;
  strconv1 needs a longer transient or a recalibrated target; examp8 needs
  a Professional reference.
- The `tsup_solve` internal test uses a folded quad4 connectivity
  (`1 2 3 4` on (0,0)(1,0)(1,1)(0,1), detJ = -0.25*eta); Gauss samples the
  fold, Lobatto did not. The corrected connectivity `1 2 4 3` passes.
  Test-input issue, not a code regression.
- The force5 fix is a general correctness fix for `force_element_edge_water`
  (independent of the quadrature): a non-selected side must not require a
  phreatic level.

## 7. Files and functions

- `polynom.cc` — `pol()`, `stress_recovery_gauss`.
- `miscel.cc` — `sri_stress_recovery_weight()`.
- `general.cc` — `general()`, the lumped sigma-dof inertia weight.
- `materi.cc` — `materi()`, the stress-dof RHS.
- `area.cc` — `area()`, the `force_element_edge_water` phreatic lookup.
- `so_bicg.cc` — `solve_iterative_bicg()`, the scale-relative CG breakdown
  test (fix A) required for dynamic8.
- `calcul_force.cc`, `delete.cc` — the replicated tensor-product default
  (`msf_element_rule`, `mc_element_rule`), switched to Gauss with the
  element default; `msf_element_rule` no longer forces the axisymmetric
  1-point rule (K20).

## Differences with the Professional version

- **`tet10` integration rule** (known difference, not decided; BACKLOG C31;
  AUDIT K04). GNU: Nodal quadrature exact only for linear integrands
  (`tet10_2` gravity sigzz -0.894). Professional 02-08-2026: Exact (-1).
  Status: Scheduled for F4 (design + corpus A/B), project owner 2026-10-06.
  Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **quad4 nodal stress recovery on distorted or axisymmetric elements**
  (known difference, not decided; BACKLOG C42). GNU: bilinear
  interpolation `1,xi,eta,xi*eta` of the 2x2 Gauss values in the natural
  coordinates (`sri_stress_recovery_weight`). Professional 02-08-2026:
  uniform-weight linear least squares `1,x,y` in the physical
  coordinates. Equal on plane rectangles; distorted quad4 3-7 %
  (investigation probe `hg4d`), axisymmetric ring `axisym2` sigzz
  10760.09 vs 10759.22 (displacements equal). A prototype (variant v4
  of F4 T4) reproduced the Professional `axisym2` value but flips the
  exit code of `cd_method`; not adopted. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`post_point_dof` of recovered fields (stress/strain)** (known
  difference, not decided; AUDIT section 5 (after the K04 table)). GNU:
  Interpolates the averaged nodal field (`validation_2` -0.036558).
  Professional 02-08-2026: Element-local recovered value from `element_dof`
  (-0.036441). Status: Both legitimate; no decision. Targets on recovered
  `post_point_dof` inherit the difference. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`number_of_integration_points`** (GNU extension / input dialect; AUDIT
  1.1). GNU: GNU-only record. Professional 02-08-2026: Unknown record; uses
  its per-operation default. Status: GNU extension; kept. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
