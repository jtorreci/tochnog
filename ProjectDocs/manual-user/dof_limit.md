# dof_limit

## Description

Input record (manual Professional 6.411):

```
dof_limit lower_dof_0 upper_dof_0 lower_dof_1 upper_dof_1 ...
```

Lower and upper allowed value of every primary dof. The pairs follow the
`dof_label` order and there must be exactly one pair per dof (all dofs,
also stresses, strains and displacements, not only the principal ones);
another number of values is an input error ("Problem reading :
dof_limit"). Use a very large range (e.g. `-1.e20 1.e20`) for a dof that
needs no limit.

After every nodal update (each equilibrium iteration and the end of the
time step) the `node_dof` value of every node is clamped to its range:
free and prescribed nodes, solved dofs and averaged dofs (stresses) alike.
A clamped velocity also limits the displacement increment
(displacement = old displacement + clamped velocity * dtime). The
integration-point values (`element_intpnt_dof`) are not clamped.

## Why it is useful

It keeps a field inside its physical range, e.g. a temperature or a
saturation that must not overshoot, or a velocity that must not exceed a
maximum.

## Parameters

- `lower_dof_i`, `upper_dof_i`: allowed range of the i-th dof of
  `dof_label` (doubles).

## Example

`validation-suite/test-2026/dof_limit1.dat`: three quad4 in a row,
`condif_temperature` only, ends prescribed at 0 and 100:

```
dof_limit 10. 40.
```

Result: node temperatures 10 33.333 40 40 (without the record 0 33.333
66.667 100). `dof_limit2.dat` limits `velx` of a rigid translation 0.01 to
0.005 (22 values: velx vely disx disy, 6 stresses, temp): node velx 0.005,
disx 0.005. Both equal the Professional 02-08-2026.

## Differences with the Professional version

- The Professional clamps the free principal dofs before the element pass
  of the step but leaves the prescribed values unclamped there (its
  element records still show 0 and 100 at the prescribed nodes), and
  clamps every node at the end. The GNU clamps every node after each
  nodal update; `bounda` re-imposes the prescribed values at the start of
  each equilibrium iteration, but the post-solve material pass of the step
  sees the clamped prescribed values. The final `node_dof` is the same; the
  element records of elements with a prescribed node outside its limit
  can differ (a prescribed value outside its own limit is a contradictory
  input). DIVERGENCES DOF-LIMIT-PRESCRIBED (class B).
- With `derivatives` the limits apply to the primary value of each dof
  only (one pair per dof, not per derivative); not compared with the
  Professional.
