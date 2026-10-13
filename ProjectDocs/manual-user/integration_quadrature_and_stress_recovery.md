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
(`-gauss`, `-lobatto`, `-reduced`, and `-minimal` / `-maximal` / `-normal`).
The reduced method and the SRI keyword are unchanged.

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
interface elements keep their nodal (Lobatto) integration. The exception is
the `tet10`: its 5-point rule would give negative corner masses, so the GNU
lumps the inertia/capacity with the positive HRZ masses (corner V/36,
mid-edge 4V/27; see the differences below).

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

- **Nodal stresses of `-lobatto` quad4 elements** (BACKLOG C79, open,
  LOW). Where elements of different volume share a node, the GNU
  averages their values at the node weighted by volume, the Professional
  takes the plain mean (two quad4 of widths 1 and 2: 6.525 vs 6.15).
  Velocities and integration-point stresses are equal. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **Lumped mass / capacity of the `tet10`** (deliberate divergence,
  project owner 2026-10-09; DIVERGENCES TET10-HRZ). GNU: positive HRZ
  lumping, corner V/36 and mid-edge 4V/27 (the diagonal of the
  consistent mass scaled to the element total). Professional 02-08-2026:
  row-sum lumping of its 5-point rule, corner -V/20 and mid-edge V/5;
  the negative corner masses make dynamics unstable (sudden gravity on a
  48-tet10 column: the Professional solver fails, a GNU row-sum build
  diverges to 5e8). Static results are unaffected. Index:
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
- **`prism6` with `group_integration_points -minimal`** (known
  difference, not decided; BACKLOG C76). The Professional uses 2 points
  (triangle centroid x 2 Gauss levels); the GNU ignores the record for a
  prism6 and keeps its 6 points.
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
