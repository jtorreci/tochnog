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

This matches the Tochnog Professional default. The `group_integration_method`
and `group_integration_points` records still override it per element group
(`-gauss`, `-lobatto`, `-reduced`, and `-minimal` / `-maximal` / `-normal`).
The reduced method and the SRI keyword are unchanged.

The mass matrix is the row-sum lumped mass (quadrature-invariant), and the
interface elements keep their nodal (Lobatto) integration.

## Nodal stress recovery

The stresses are stored at the nodes. For the interior Gauss rules the
GNU version now recovers them the way the Professional does: a per-element
least-squares fit of the integration-point stresses evaluated at the nodes,
assembled as the arithmetic mean over the elements that share the node. The
resulting nodal stress field matches the Professional (the `beam2d` /
`beam2d_1` cantilever corner stress and the `post4` moment are the
discriminating examples).

## When a stored reference value changes

Because the default quadrature changed, stored static/dynamic reference
values produced under the old Lobatto default (for example some test
targets calibrated to it) can change. Re-run the model and recalibrate the
target if the new value matches the Professional.

## Differences with the Professional version

- **`tet10` integration rule** (known difference, not decided; BACKLOG C31;
  AUDIT K04). GNU: Nodal quadrature exact only for linear integrands
  (`tet10_2` gravity sigzz -0.894). Professional 02-08-2026: Exact (-1).
  Status: Scheduled for F4 (design + corpus A/B), project owner 2026-10-06.
  Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **axisymmetric hoop (zz) values, `materi_velocity` quad4 (`axisym2`,
  `axisym9`)** (known difference, not decided; AUDIT K20; BACKLOG C32). GNU:
  1-point Gauss for axisymmetric `materi_velocity` quad4: element-constant
  hoop term (axisym2 1.000e-3; axisym9 0.00158902). Professional 02-08-2026:
  Linear hoop field (1.0498e-3; axisym9 0.0015915). Status: Scheduled for F4
  (design + corpus A/B), project owner 2026-10-06. Index:
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
