# element_dof

## Description

Output record (manual Professional 6.421): "Unknowns as saved per element in
the element nodes. First dof's in the first node. Then dof's in the second
node. Etc." One record per element:

```
element_dof <element> u_1(n1) ... u_n(n1) u_1(n2) ... u_n(n2) ...
```

`nnol*n` values (`nnol` nodes of the element, `n` unknowns in `dof_label`
order), filled after every converged time step for every element with an
integration-point record (`element_intpnt_dof`). Per node:

- principal unknowns (velocities, temperature, pressure, ...): the value of
  that node (`node_dof`);
- the other unknowns (strains, stresses, kappa, damage, history variables):
  the nodal value recovered from the integration points of THIS element
  only (the per-element fit of the Professional: the point value for the
  1-point rules and for Lobatto rules with the points at the nodes, a
  linear least-squares fit in physical coordinates otherwise, e.g. the
  quad4 2x2 Gauss rule).

`node_dof` averages these per-element values over the elements of a node,
so on shared nodes `element_dof` and `node_dof` differ.

Output only: it is never read by the calculation, and it is not meant to be
given in the input file.

## Example

`element_dof_nodal1.dat` (validation-suite/test-2026), two quad4: at the
shared node 2 element 1 gives sigxx -1.6178, element 2 -3.0275, and
`node_dof` -2.3227 (the mean). Also `element_dof_nodal_1pt1.dat` (1-point
rule: the point value at every node) and `element_dof_nodal_quad9_1.dat`
(quad9 with 2x2 points). All equal the Professional 02-08-2026 except the
displacement slots.

## Differences with the Professional version

- The displacement slots (`disx`, `disy`, `disz`) hold the recovery of the
  GNU integration-point displacement; the Professional writes 0 there
  (DIVERGENCES ELEMENT-DOF-DISPLACEMENT, class B).
- Rules without a measured Professional fit use a GNU fallback: the mean of
  the points weighted with the shape function of the node and the point
  volume. Every rule `pol()` currently produces has a fit; the fallback is a
  guard.
- Format change of the GNU .dbs (F5 T2 V5, 2026-11-02): until then the GNU
  `element_dof` held the integration-point values (now
  [element_intpnt_dof](element_intpnt_dof.md)), padded to 27 points. Scripts
  that read the old layout must read `element_intpnt_dof` now.
