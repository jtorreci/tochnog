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
`tet10` uses the Professional linear least-squares fit `1,L1,L2,L3` of its
5 points. The resulting nodal stress field matches the Professional (the `beam2d` /
`beam2d_1` cantilever corner stress and the `post4` moment are the
discriminating examples).

The recovery factors of quad4, quad9, hex8, hex27 and tet10 are computed
once per run (F4 T9, 2026-10-10): the results are unchanged and large 3D
models recover their nodal stresses faster (Professional `large1`
37 s -> 27 s).

## When a stored reference value changes

Because the default quadrature changed, stored static/dynamic reference
values produced under the old Lobatto default (for example some test
targets calibrated to it) can change. Re-run the model and recalibrate the
target if the new value matches the Professional.

## Differences with the Professional version

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
  Professional. Still different: the 3x3 quad9 and the hex8 on curved or
  distorted elements (GNU natural-coordinate fit, BACKLOG C70), and the
  convected-stress models (Eulerian stress path not consistent with the
  Professional, BACKLOG C72; same recovery operator since 2026-10-12).
  Index:
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
