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
