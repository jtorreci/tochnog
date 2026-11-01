# Integration quadrature and nodal stress recovery

This note describes the default element integration of the GNU version and
how the nodal stresses are recovered from the integration-point stresses.
It exists because the GNU version aligns these with Tochnog Professional.

## Default integration

For the tensor-product elements (`bar2`, `bar3`, `bar4`, `quad4`, `quad9`,
`quad16`, `hex8`, `hex27`, `hex64`) the default is the **Gauss** rule with
one point per node direction:

| element | default rule |
|---|---|
| quad4 / hex8 | 2x2 / 2x2x2 Gauss |
| quad9 / hex27 | 3x3 / 3x3x3 Gauss |
| tet10 | Keast 5-point rule (centroid + 4 interior points, degree 3) |

This matches the Tochnog Professional default. Axisymmetric element groups
(`group_axisymmetric -yes`) use the same default since F4 T4 (2026-10-09):
the quad4 is integrated with 2x2 Gauss and the quad9 with 3x3 Gauss, as the
Professional does. Before, the GNU forced a 1-point rule for axisymmetric
`materi_velocity` and refused the quadratic elements ("not available for
axisymmetric analysis"); the hoop term was then constant over the element
and the quad4 hourglass modes were not restrained.

The `group_integration_method`
and `group_integration_points` records still override it per element group
(`-gauss`, `-lobatto`, `-reduced`, and `-minimal` / `-maximal` /
`-normal`). The `-reduced` method with its factor record is described
in [group_integration_method_reduced](group_integration_method_reduced.md)
(F4 T6f C80, 2026-10-16).
The SRI keyword is unchanged.

With `group_integration_method -lobatto` in an axisymmetric group, an
element that touches the symmetry axis (a node with radius below 1e-9)
is integrated with the Gauss rule instead, as the Professional does:
the Lobatto points at r=0 have zero volume and an infinite hoop strain
term u/r (before 2026-10-13 the GNU produced NaN velocities, sfnet
`tube1`). The other elements of the group keep the Lobatto rule.
Minimal example (one element on the axis):

```
group_type 0 -materi
group_axisymmetric 0 -yes
group_integration_method 0 -lobatto
node 1 0. 0.
node 2 1. 0.
node 3 0. 1.
node 4 1. 1.
element 1 -quad4 1 2 3 4
```

The `tet10` rule is the Professional one since F4 T4 (2026-10-09). The
former nodal rule integrated only linear integrands exactly, so a
quadratic tetrahedron mesh under gravity gave sigzz -0.894 instead of the
exact -1 (`tet10_2`); the 5-point rule gives -1.

The mass matrix is the row-sum lumped mass (quadrature-invariant), and the
interface elements keep their nodal (Lobatto) integration. The exceptions
are the `tet10` and the `tria6`: their row sums give negative (tet10) or
zero / negative (tria6) corner masses, so the GNU lumps them with the
positive HRZ masses (the consistent-mass diagonal scaled to the element
total: tet10 corner V/36, mid-edge 4V/27; tria6 corner A/19, mid-side
16A/57, r-weighted when axisymmetric; see the differences below).

The `tria6` (since 2026-10-16, F4 T6f C66) integrates with the manual
6.623 rules, = Professional: 7 points (degree 5) by default and with
`-maximal`, 4 points (degree 3, centroid weight -27/48) with `-minimal`;
axisymmetric groups are accepted, and an element with a node on the axis
keeps its rule (all points are interior). Before, a 6-point nodal rule was
used for every value and axisymmetric tria6 groups were refused.

**Limitation of the tria6 `-minimal` rule** (owner decision 2026-10-17):
the negative centroid weight makes the element stiffness indefinite (one
negative eigenvalue) on curved or axisymmetric elements once Poisson's
ratio exceeds about 0.39-0.49 (plane strain curved 0.43-0.49, axisymmetric
on the axis 0.39). Straight plane elements are exact (equal to the
7-point stiffness) and the default 7-point rule is positive semi-definite
in every case. Do not use `-minimal` for tria6 with curved or
axisymmetric elements and nu > 0.4; keep the default rule. The prism15
`-minimal` 8-point rule of the Professional is refused instead (C82)
because it is indefinite already on a regular prism at nu = 0.3.

## Nodal stress recovery

The stresses are stored at the nodes. For the interior Gauss rules the
GNU version now recovers them the way the Professional does: a per-element
least-squares fit of the integration-point stresses evaluated at the nodes,
assembled as the arithmetic mean over the elements that share the node. The
fit is the linear field `1,x,y` (2D) or `1,x,y,z` (3D) in the PHYSICAL
coordinates for every quad/brick Gauss rule (quad4 2x2, quad9 2x2 and
3x3, hex8 2x2x2, hex27 2x2x2 and 3x3x3), the `tet10` 5-point rule and
the `prism6` 6-point rule; the 1-point quad4/hex8 take the point value
(F4 T6f K21 2026-10-12 and C70 2026-10-13). The resulting nodal stress
field matches the Professional (the `beam2d` / `beam2d_1` cantilever
corner stress and the `post4` moment are the discriminating examples;
`test-2026/recovery_*` cover distorted and curved elements).

With `group_integration_method -lobatto` (F4 T6f C79, 2026-10-13) the
same arithmetic mean over the elements applies. When the Lobatto points
sit at the nodes (quad4, quad9, hex8, hex27 with the default
`-maximal` points) each element contributes the value of its point at
the node; with `group_integration_points -minimal` the quad9 / hex27
points are the element corners and each element contributes the
physical linear fit of those corner values, as for its 2x2 Gauss rule.
Before, the GNU weighted the Lobatto elements by their volume (a node
shared by a small and a large element leaned towards the large one) and
gave zero stress at the mid-side and centre nodes of a `-minimal`
quad9 / hex27. Minimal example (two elements of different size sharing
the nodes 2 and 5):

```
group_type 0 -materi
group_integration_method 0 -lobatto
node 1 0. 0.
node 2 1. 0.
node 3 3. 0.
node 4 0. 1.
node 5 1. 1.
node 6 3. 1.
element 1 -quad4 1 2 4 5
element 2 -quad4 2 3 5 6
```

The `post_calcul -materi_stress -force` section forces of an
axisymmetric `-lobatto` element touching the axis use its Gauss rule
too (the rule it was integrated with; before, the Lobatto rule gave a
wrong moment).

The recovery operator depends on the element geometry, so it is built
per element and cached while the element is visited (since C70 also for
the 3D elements; before, the natural-coordinate factors were computed once
per run, F4 T9).

## When a stored reference value changes

Because the default quadrature changed, stored static/dynamic reference
values produced under the old Lobatto default (for example some test
targets calibrated to it) can change. Re-run the model and recalibrate the
target if the new value matches the Professional.

## Differences with the Professional version

- **Nodal stresses of `-lobatto` elements** (BACKLOG C79) — RESOLVED
  2026-10-13 (F4 T6f C79): quad4, quad9, hex8 and hex27 Lobatto rules now
  = Professional (arithmetic mean; sfnet `tube1` node (1.25,0) sigxx
  -407.64, before -494.42). The bar2/bar3/bar4 and quad16 rules (Gauss,
  Lobatto, `-reduced`) follow the Professional since 2026-10-16 (BACKLOG
  C67, [nodal_stress_recovery](nodal_stress_recovery.md)); the `hex64`
  is accepted by neither version.
- **`group_integration_method -reduced`** (BACKLOG C80) — implemented
  2026-10-16 (F4 T6f C80), = Professional when the factor is given.
  Without `group_integration_method_reduced_factor` the GNU uses the
  manual default 10 and the Professional uses 1 (deliberate divergence,
  class A, project owner 2026-10-16). Details and the hourglass warning:
  [group_integration_method_reduced](group_integration_method_reduced.md).
- **Lumped mass / capacity of the `tet10`** (deliberate divergence,
  project owner 2026-10-09; DIVERGENCES TET10-HRZ). GNU: positive HRZ
  lumping, corner V/36 and mid-edge 4V/27 (the diagonal of the
  consistent mass scaled to the element total). Professional 02-08-2026:
  row-sum lumping of its 5-point rule, corner -V/20 and mid-edge V/5;
  the negative corner masses make dynamics unstable (sudden gravity on a
  48-tet10 column: the Professional solver fails, a GNU row-sum build
  diverges to 5e8). Static results are unaffected.
- **Lumped mass / capacity of the `tria6`** (deliberate divergence,
  project owner 2026-10-16; DIVERGENCES TRIA6-HRZ). GNU: positive HRZ
  lumping of every h-weighted dof (velocity, temperature, pressure,
  density, the h-weighted state variables): plane corner A/19 and
  mid-side 16A/57; axisymmetric the r-weighted diagonal int N_i^2 2 pi r
  dA scaled to 2 pi A r_c; the same masses for every
  `group_integration_points` value. Professional 02-08-2026: row sum,
  corner 0 in plane analysis (a free one-element dynamic step fails with
  a zero pivot), negative corners in axisymmetric analysis (-pi/60 at
  r = 1 for r = 1..2). Static results are unaffected. The groundflow
  velocity and the condif absorption of a `tria6` use the same HRZ
  weight since C85 (2026-10-17): a uniform pressure gradient gives the
  velocity k dp/dx at every node, corners included, and a uniform
  temperature decays at the same rate at every node. Tests
  `tria6_mass_hrz1`, `tria6_mass_hrz_axisym1`,
  `tria6_groundflow_velocity_hrz1`, `tria6_condif_absorption_hrz1`,
  `tria6_condif_absorption_hrz_axisym1`. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`post_point_dof` inside a `tet10`** (deliberate divergence, project
  owner 2026-10-09; DIVERGENCES TET10-PPOINT). GNU: interpolates the
  nodal values with the quadratic `tet10` shape functions (all 10
  nodes), the consistent finite-element interpolation. Professional
  02-08-2026: interpolates linearly from the 4 corner nodes (`tet10p`
  post points differ by up to 2 %). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **quad4 nodal stress recovery on distorted or axisymmetric elements**
  (BACKLOG C42) — RESOLVED 2026-10-12 (F4 T6f K21): the 2x2 Gauss quad4
  (and the 2x2 quad9) now use the Professional linear least squares
  `1,x,y` in the physical coordinates (on-axis axisymmetric elements:
  shape-function weighted average); `axisym2` inner sigzz 10759.2155 =
  Professional. The 3x3 quad9, hex8, hex27, tet10 and prism6 followed on
  2026-10-13 (F4 T6f C70; before: natural-coordinate fits, different on
  curved or distorted elements). Still different: the prism6
  integration points (deliberate, see the next item) and the convected-stress models (BACKLOG C72: since 2026-10-13 the nodal
  stress is recovered from the integration points as in the
  Professional; Professional-side residuals remain on `validation_1`).
  Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`prism6` integration points** (deliberate divergence, project owner
  2026-10-13, BACKLOG C75). The Professional integrates a prism6 at the
  three triangle CORNERS x 2 Gauss levels in zeta (weight 1/12 each; also
  with `-maximal` and with inertia); the GNU keeps the interior points
  (2/3,1/6,1/6) x 2 Gauss levels, exact to degree 2 in the triangle. The
  GNU rule halves the bending error of a prism6 cantilever (40x2x4 mesh:
  -4.85 % against -9.95 % of the analytic deflection). Results agree when
  the stress is linear in the element; otherwise prism6 stresses and
  displacements differ slightly (distorted-prism probe: nodal sigzz 1.6 %
  of the scale).
- **`prism6` with `group_integration_points -minimal`**: = Professional
  since 2026-10-16 (C76, closed by following it): 2 points, the triangle
  centroid x 2 Gauss levels, weight 1/4 each. An under-integrated rule
  (1 spurious mode per element, restrained by the mesh in every case
  tried) that halves the bending error of the 6-point rule (cantilever
  20x2x2: -7.1 % against -13.1 %). Before, the GNU ignored the record.
- **`prism15`**: the 21-point rule had two wrong triangle points (one
  outside the triangle) until 2026-10-16; it is now the Professional
  degree-5 rule (elastic displacements 3e-11 from the Professional). With
  `-minimal` the Professional uses 8 points (a 4-point triangle rule with a
  negative weight) whose element stiffness is indefinite; the GNU keeps 21
  points (deliberate divergence, owner decision 2026-10-16; BACKLOG C82).
- **`tet4` with `-maximal`**: 5 points (the Keast rule of the tet10), as
  manual 6.623 and the Professional, since 2026-10-16 (before: 4 points
  at the nodes; same stiffness for a linear material).
- **`tria6` rule**: = manual 6.623 and Professional since 2026-10-16
  (F4 T6f C66): 7 points (default, `-maximal`), 4 points (`-minimal`),
  axisymmetric accepted; nodal stresses equal within 3e-10 of the scale.
  Its lumped mass is the class A row below.
- **`tria3` `-maximal`, `tet10` `-maximal`** (deliberate divergence from
  the MANUAL, owner decision 2026-10-16; BACKLOG C83): the manual (6.623)
  specifies 4 and 10 points but gives neither points nor weights; both
  versions use 3 (the nodes) and 5 (Keast), so the GNU equals the
  Professional. The classic 4-point triangle rule has the negative
  centroid weight rejected for the prism15 (C82); a positive 10-point tet
  rule exists (Keast) but no version uses it.
- **Axisymmetric `tria3` `-maximal`** (= Professional since 2026-10-16;
  BACKLOG C83): 3 points at the nodes; an element with a node on the axis
  uses 1 point at the centroid (the hoop term of a point at r = 0 is
  infinite). Before, the GNU used 1 point for every axisymmetric tria3.
  Default and `-minimal`: 1 point in both versions.
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
