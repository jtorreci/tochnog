# element_beam_force_moment

## Description

`element_beam_force_moment <element> force_x_first_node force_y_first_node
force_z_first_node moment_x_first_node moment_y_first_node
moment_z_first_node force_x_second_node ... moment_z_second_node` (manual
Professional 6.416): after the calculation this record holds the forces
and moments of a beam in the local beam axes x, y, z (12 values: first
node, then second node).

The values are the end forces that the nodes apply on the beam, so the
values at the first node have a minus compared with the second node: a
beam with a constant z moment of 10 gives

```
element_beam_force_moment 20 0 0 0 0 0 -10 0 0 0 0 0 10
```

It is an output record, written for every `-beam` and `-trussbeam`
element; you do not give it in the input.

- Beam with local axes (a `group_beam_direction_z`,
  `group_beam_direction_z_reference_point` or `element_beam_direction_z`
  record applies, see [group_beam_direction_z](group_beam_direction_z.md)):
  all 12 components in the local axes.
- Classic GNU beam (no such record): the beam is planar, so only
  force_x, force_y and moment_z are filled, in the axes of its plane
  (local x along the beam, local y in the plane, local z out of the plane,
  the global z axis for the default x-y plane). The other components are 0.

The beam has no axial stiffness: force_x is 0 (up to round-off); the
axial force of a `-trussbeam` is in `element_truss_force`.

`control_print_beam_force_moment` uses this record for the beams with
local axes.

## Example

`element_beam_force_moment1.dat` (validation-suite/test-2026): a 3D
cantilever of length 1 along x, `group_beam_young 0 1`,
`group_beam_inertia 0 1 4 1`, tip moment 10 about z (`bounda_force 1 2
-rotz`, `bounda_time 1 10.`): `element_beam_force_moment 1 0 0 0 0 0 -10
0 0 0 0 0 10`, tip vely 1.25, rotz 2.5 (the manual example).

`element_beam_force_moment2.dat`: out-of-plane force, torsion and a
moment about y on beams with local axes.

## Differences with the Professional version

- A beam whose axial displacement is free stops the GNU direct solver;
  the Professional solves it (BEAM-AXIAL-FREE). Restrain the axial dof or
  add a truss.
- Default `group_beam_memory` (GNU `-updated`, Professional
  `-total_linear`) changes the values slightly for finite displacements
  (BEAM-MEMORY-DEFAULT).
- 3D beam without a local z record: only the in-plane components (C37).
- The Professional 2D beam could not be measured (singular on the test
  inputs); the GNU 2D values are checked against the analytic answer
  (`element_beam_force_moment_2d1`).
