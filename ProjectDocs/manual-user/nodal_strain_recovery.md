# nodal_strain_recovery — nodal values of the strain tensors

## What it does

The strain dofs `materi_strain_total` (`ept..`), `materi_strain_elasti`
(`epe..`) and `materi_strain_plasti` (`epp..`) are integration-point
quantities written to the nodes. Since 2026-10-06 (audit K11) they are
recovered with the SAME operator as the stress: the Professional
per-element least-squares fit of the integration-point values (bilinear
interpolation of the 2x2 Gauss values for quad4, linear fit for
quad9/hex8/hex27; the plain shape-function weight for bar/tria/tet and
the Lobatto rules), the nodal value being the mean of the adjacent
elements. Before, the strain used a shape-function-weighted average while
the stress already used the Professional operator, so strain and stress
nodal fields were not consistent (and both differed at boundaries from
the Professional).

## Why it is useful

Nodal strains now equal the Professional and are consistent with the
nodal stresses: on a clamped edge the recovered `eptyy` is exactly 0
(the Professional value), where the old average gave -8.6e-5 in the test
below; the shear strain near the clamp is no longer smeared (factor 3 in
`tslv_part`). The element/integration-point values and the solution are
unchanged.

## Input

Nothing new; the recovery applies whenever the strain dofs are declared
in the initia part. Elements with selective reduced integration
(`group_materi_sri`) keep the shape-function weight.

## Example

`validation-suite/test-2026/strain_recovery_clamped1.dat`: two quad4,
E 1, nu 0.3, left edge clamped, right edge pulled 1e-3:
eptyy(node 1) 0, eptxy(node 1) 6.44824708e-5, eptyy(node 2)
-2.57929883e-4 (= Professional 02-08-2026).

## Differences with the Professional version

None measured on the corpus inputs the Professional runs (14 tests move
to the Professional values; the other changes are round-off of zero
fields).
