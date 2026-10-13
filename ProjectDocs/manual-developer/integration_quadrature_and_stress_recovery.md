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
| element stiffness, tet10 | **Keast 5-point** (centroid -4/5, (1/2,1/6,1/6,1/6) x4 9/20) | nodal rule (degree 1) | Keast 5-point (F4 T4) | `tet10_2` gravity sigzz -1; `.dbs` npoint 5 |
| mass matrix | **row-sum lumped** (quadrature-invariant) | row-sum lumped | unchanged; tet10 HRZ (F4 T4, deliberate) | distorted-element dynamic probe; `mass10` |
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
| quad4 | 4 | ~~bilinear interpolation `1,x,y,xy`~~ -> **linear least squares `1,x,y` in PHYSICAL coordinates** (K21, 2026-10-12; axisymmetric with a node at r=0: `sum h sig / sum h`) | 1e-12 (distorted probes, `tube1`) |
| quad4 | 1 | the point value (K21) | exact (`examp1`) |
| quad9 | 4 | linear least squares `1,x,y` in physical coordinates (K21) | 1e-12 (`incnav4`, probes) |
| quad9 | 9 | **uniform-weight linear least squares** `1,x,y` (3 terms) in PHYSICAL coordinates (C70, 2026-10-13; before: GNU table in natural coordinates); axisymmetric with a node at r=0 (9 and 4 points): `sum h sig / sum h` | 1e-12 (curved probes) |
| hex8 | 8 | **uniform-weight linear least squares** `1,x,y,z` (4 terms) in PHYSICAL coordinates (C70; before: natural table) | 1e-12 (distorted probes) |
| hex8 | 1 | the point value (C70; before: volume-weighted h) | exact |
| hex27 | 27 / 8 | uniform-weight linear least squares `1,x,y,z` in PHYSICAL coordinates (C70, measured; before: natural table for 27, h for 8) | 1e-12 (curved probes) |
| tet10 | 5 | uniform-weight linear least squares `1,x,y,z` in PHYSICAL coordinates (C70; the F4 T4 barycentric fit `1,L1,L2,L3` is the same on straight-edged tet10 only) | 1e-12 (curved probes); `tet10_2` exact (-1) |
| prism6 | 6 | uniform-weight linear least squares `1,x,y,z` in PHYSICAL coordinates (C70; before: volume-weighted h). The Professional points are the triangle corners x 2 Gauss levels, the GNU points (2/3,1/6,1/6) x 2 Gauss levels (C75) | 1e-13 |
| tria3 / tet4 / bar2 1 point, tria3 3-point and tria6 6-point nodal rules | 1 / 3 / 6 | the point value / the value of the point at the node (F4 T6f, two-element probes) | exact |

The operator was believed FIXED in the element natural coordinates (the
same 9x9 matrix fits all 45 beam2d_1 elements of different geometry to
1e-9: they are parallelograms, where the physical and natural fits
coincide). K21 (2026-10-12) showed it is a fit in the PHYSICAL
coordinates of the points (distorted quad4, curved quad9, distorted hex8
probes); the GNU builds it per element for every Gauss rule since C70
(2026-10-13; K21 did the 2D 2x2 rules first, `nodal_stress_recovery.md`). The
least-squares fit uses UNIFORM weights (the Gauss weights give a different
operator). The order rule is not a single clean formula: quad4 uses the
full 4-term interpolation while quad9 and hex8 use the linear fit. The
operators above were determined empirically and are the discriminating
data; hex27/quad16 were not independently exercised.

Resolved 2026-10-12 (K21; see the table above and
`nodal_stress_recovery.md`). Historical note, the open question of
F4 T6f (2026-10-11, BACKLOG C64):
on a Drucker-Prager variant of `validation_14` the Professional
per-element quad4 values (`element_dof` from `element_intpnt_dof`) are
reproduced exactly in 92 of 95 elements by a LINEAR least-squares fit
`1,xi,eta` of the 2x2 Gauss values, not by the bilinear interpolation
`1,x,y,xy` of the table above. The two operators agree whenever the
Gauss values have no `xi*eta` mode, so the finding discriminates only if
those elements carry that mode; the original identification used probes
whose stress field may not have had it either. An elastic quad4 probe with a non-zero `xi*eta` stress mode is
needed before changing the operator (see also C42, physical vs natural
coordinates).

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

### 3.1b Simplex elements and the 1-point bar2 (F4 T6f, 2026-10-11)

`sri_stress_recovery_weight()` returned `h_inol` for the tria3, tria6,
tet4 and bar2, so their nodes were the area/volume/length-weighted
average of the adjacent elements. The Professional 02-08-2026 uses the
arithmetic mean for them too (two-element probes, exact). The weight is
now `F/volume` (`F = 1` for the 1-point rules, `F = h` for the tria3
3-point and tria6 6-point nodal rules), like the table elements. Unmeasured
rules (tet4 `-maximal`, bar2 2 points, bar3/bar4, the quad Lobatto
rules) keep `h_inol` (BACKLOG C67). Details and tests:
`nodal_stress_recovery.md`.

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
constant over the element. The remaining `axisym2` sigzz difference was
the quad4 nodal recovery (natural-coordinate bilinear interpolation in
the GNU, physical-coordinate linear least squares in the Professional;
BACKLOG C42, resolved by K21 2026-10-12: `axisym_quad4_ring1` sigzz
10759.2155 = Professional). `msf_axisym` (`post_calcul -materi_stress -force`): the
section moment changes from 0 to 1/960, the analytic moment of the
2*PI*r-weighted face (the build_safe tolerance 0.005 holds).
`foundation_consolidation_settlement` (Pro `other`, axisymmetric) now
runs past the point where the 1-point rule failed and needs more than
4 GB (F4 T7 memory hogs).

### 3.3 tet10 rule, recovery and lumping (C31, F4 T4, 2026-10-09)

- `polynom.cc` `pol()`, TET10 branch: `npoint = 5`, Keast rule (centroid
  weight -4/5, the four points with one barycentric coordinate 1/2 and
  the others 1/6, weight 9/20; normalised to 1 like every rule of
  `pol()`), `stress_recovery_gauss = 1`. Degree 3: the quadratic
  stiffness (degree 2) and the consistent body load (corner -V/20,
  mid-edge V/5) are exact. The former nodal rule (weight EPS_TRIANGLE
  at the corners, ~1/6 at the mid-edge nodes) was exact only for linear
  integrands.
- `miscel.cc` `sri_stress_recovery_weight()`: tet10 branch. The
  Professional fits `1,L1,L2,L3` to the 5 point values with uniform
  weights. In the barycentric basis the normal matrix is
  `I/9 + (41/144) J`, so the factor is the closed form
  `F = 9 L(node).L(point) - 41/20` (sums to 1 over the points); the
  volume cancellation of section 3 applies (the centroid point volume is
  negative, the cancellation still gives `volume*weight = F`).
- `general.cc` `general()`, lumped inertia: for a 5-point tet10 and a
  principal dof (`principal_unknown`: velocity, temperature, groundflow
  pressure, density, wave) the weight is the HRZ mass fraction, corner
  1/36 and mid-edge 4/27 (consistent-mass diagonal 6/420 and 32/420
  scaled to the total; sum 1). The integration-point volumes sum to V,
  so the node gets V/36 or 4V/27. The non-principal recovered dofs keep
  `h`: their `materi.cc` numerators are `h`-weighted, so numerator and
  denominator stay consistent (corner weights are all negative and
  mid-edge all positive, both proper averages). The stress/strain dofs
  keep the recovery weight. HRZ is a DELIBERATE DIVERGENCE (owner
  2026-10-09): the Professional row sum gives the corner -V/20.
- Probes (Professional 02-08-2026 vs GNU, F4 T4): `tet10_2` sigzz
  -0.894 -> -1.00000000515 (Pro -1); `tet10_1` stays exact; `mass10`
  (one tet10, unit point forces, dt 0.01): GNU velocities 2.16 corner /
  0.405 mid-edge = dt/m with the HRZ masses (Pro -1.2 / 0.3, row sum);
  `dyn10` (48 tet10 column, sudden gravity, inertia): GNU stable, top
  disz -0.5286 at t 3 (analytic undamped -0.5), Pro solver failure,
  GNU 5-point build with row sum 5e8.
- Tests: `test-2026/tet10_body_load1`, `tet10_mass_hrz1`,
  `tet10_dynamic_gravity1`.

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
- quad4 operator: resolved by K21 (2026-10-12, C64 closed): the
  Professional fits `1,x,y` in physical coordinates for every measured
  type; the 92/95 `validation_14` match was the natural fit agreeing with
  the physical one on near-parallelograms. The 3x3 quad9 and the 3D
  types followed in C70 (2026-10-13, all measured physical). Unmeasured:
  prism15, quad16/hex64, the Lobatto rules (C67).
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

## 6b. Recovery factor table (F4 T9, 2026-10-10)

The least-squares factor `F(inol, ipoint)` depends only on the element
type, the node and the integration point, but
`sri_stress_recovery_weight()` rebuilt and inverted the 4x4 normal matrix
on every call - once per node, integration point and element on every
recovery and inertia pass (43.7 million fits on the Professional
`large1`). `miscel.cc` now builds a static table once
(`sri_recovery_table_build()`, a function-local `static const`, so the
first call initialises it thread-safely) for the five recovered types:
quad4 (4 points, bilinear), quad9 (9, linear), hex8 (8, linear), hex27
(27, linear) and tet10 (5 Keast points, closed form). The fit itself
moved unchanged into `sri_recovery_factor()` (the space dimension is
the type's own instead of the global `ndim`, which the type test already
fixed); a singular normal matrix keeps the `h_inol` fallback through the
`ok` flag. The SRI Lagrange path and the volume cancellation are
unchanged.

Results are bit-identical: every suite `.dbs` equal to the old build
(297 files), corpus A/B on omen (F4 T9 SEGUIMIENTO row). Time: `large1`
37 s -> 27 s locally (4 GB cap).

Superseded by C70 (F4 T6f, 2026-10-13): the Professional operator depends
on the element geometry, so the table, `sri_recovery_factor()` and
`sri_recovery_table_build()` were removed; `recovery_fit_factor()` builds
the operator once per element visit from the node coordinates and caches
it per thread (`nodal_stress_recovery.md`), which keeps the cost at one
small fit per element instead of one per (node, point) call.

## 7. Files and functions

- `polynom.cc` — `pol()` (incl. the tet10 Keast rule, C31),
  `stress_recovery_gauss`.
- `miscel.cc` — `sri_stress_recovery_weight()`, `recovery_fit_kind()`,
  `recovery_fit_shape()`, `recovery_fit_build()`, `recovery_fit_factor()`
  (K21 / C70 per-element physical fit; the F4 T9 factor table was
  removed in C70).
- `general.cc` — `general()`, the lumped sigma-dof inertia weight and the
  tet10 HRZ lumping (C31).
- `materi.cc` — `materi()`, the stress-dof RHS.
- `area.cc` — `area()`, the `force_element_edge_water` phreatic lookup.
- `so_bicg.cc` — `solve_iterative_bicg()`, the scale-relative CG breakdown
  test (fix A) required for dynamic8.
- `calcul_force.cc`, `delete.cc` — the replicated tensor-product default
  (`msf_element_rule`, `mc_element_rule`), switched to Gauss with the
  element default; `msf_element_rule` no longer forces the axisymmetric
  1-point rule (K20).

## Differences with the Professional version

- **Lumped mass / capacity of the `tet10`** (deliberate divergence,
  project owner 2026-10-09; DIVERGENCES TET10-HRZ). GNU: positive HRZ
  lumping, corner V/36 and mid-edge 4V/27 (the diagonal of the
  consistent mass scaled to the element total). Professional 02-08-2026:
  row-sum lumping of its 5-point rule, corner -V/20 and mid-edge V/5;
  the negative corner masses make dynamics unstable (sudden gravity on a
  48-tet10 column: the Professional solver fails, a GNU row-sum build
  diverges to 5e8). Static results are unaffected. Code: `general.cc`
  `general()` (HRZ weight). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`post_point_dof` inside a `tet10`** (deliberate divergence, project
  owner 2026-10-09; DIVERGENCES TET10-PPOINT). GNU: interpolates the
  nodal values with the quadratic `tet10` shape functions (all 10
  nodes), the consistent finite-element interpolation. Professional
  02-08-2026: interpolates linearly from the 4 corner nodes (`tet10p`
  post points differ by up to 2 %). Code: `point_el.cc` (tet10 shape functions). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **quad4 nodal stress recovery on distorted or axisymmetric elements**
  (BACKLOG C42) — RESOLVED 2026-10-12 (F4 T6f K21): `miscel.cc`
  `recovery_fit_factor()` builds the Professional physical-coordinate
  linear fit per element for the 2x2 quad4/quad9 (on-axis `h` rule,
  1-point quad4 point value); `axisym_quad4_ring1` sigzz 10759.2155 =
  Professional. The convected-stress gate of WU1 was removed in WU1b
  (2026-10-12): the `cd_method` flip came from an ill-posed input
  (`-ra 1 4 -ra` is a node list), rewritten as a well-posed shear; the
  convected stress path is BACKLOG C72 (nodal output recovered
  absolutely since F4 T6f C72, `c1a1ae9`; Professional-side residuals
  class B). The 3x3 quad9, hex8, hex27, tet10 and prism6 followed in C70
  (2026-10-13). The prism6 integration points stay different on purpose
  (see the next item). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`prism6` integration points** (deliberate divergence, project owner
  2026-10-13, BACKLOG C75, DIVERGENCES class A). Professional 02-08-2026:
  triangle corners (1,0),(0,1),(0,0) x zeta = 1/2 +- 1/(2 sqrt 3),
  weight 1/12 (`element_intpnt_iso_coord`; default, `-normal`,
  `-maximal`, `inertia_apply -yes`). GNU `pol()` (`polynom.cc`, PRISM6
  branch): (2/3,1/6,1/6) x the same levels, weight 1/12. Both full rank
  and patch-test exact; the vertex rule integrates L_i^2 as A/3 instead
  of A/6 and stiffens bending (cantilever 40x2x4: -9.95 % against
  -4.85 %). `recovery_fit_shape` (`miscel.cc`) repeats the GNU point
  table and must stay consistent with `pol()`.
- **`prism6` `-minimal`** (BACKLOG C76, class B): Professional 2 centroid
  points; the GNU `pol()` reads `GROUP_INTEGRATION_POINTS` but does not
  branch on it for a prism6, and `recovery_fit_kind` accepts only
  npoint 6. Not implemented.
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
