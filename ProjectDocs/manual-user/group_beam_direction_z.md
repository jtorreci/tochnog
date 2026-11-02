# group_beam_direction_z + group_beam_direction_z_reference_point + element_beam_direction_z

## Description

These records define the local beam axes (manual Professional 6.593,
6.594 and 6.415):

- `group_beam_direction_z <group> dir_z,x dir_z,y dir_z,z`: the local
  beam z direction, in global axes, for every beam of the element group.
  The vector must be perpendicular to the beam length direction (it does
  not need to be a unit vector).
- `group_beam_direction_z_reference_point <group> point_x point_y
  point_z`: the local z axis points to this point. The vector from the
  beam middle point to the point is taken; its part perpendicular to the
  beam is the local z axis. Useful for a tunnel lining: give the middle of
  the tunnel axis and every lining beam gets its z axis towards it.
- `element_beam_direction_z <element> dir_z,x dir_z,y dir_z,z`: the same
  as `group_beam_direction_z`, for one beam element.

Priority: `element_beam_direction_z`, then `group_beam_direction_z`, then
`group_beam_direction_z_reference_point`.

The local axes are: x = from the first to the second node of the beam,
z = the direction above, y = z cross x (right handed).

When one of the records applies to a beam, the beam is a full 3D beam in
its local axes:

- bending in the local x-y plane with `E * Izz`,
- bending in the local x-z plane with `E * Iyy`,
- torsion about the local x axis with `G * J`,

with `E` from `group_beam_young`, `G` from `group_beam_shear` and
`Iyy Izz J` from `group_beam_inertia`. As in the Professional the beam
has no axial stiffness: combine it with a truss (or a `-trussbeam`
element) for the axial force, or restrain the axial displacement. The
forces and moments are written in `element_beam_force_moment` (see
[element_beam_force_moment](element_beam_force_moment.md)).

The records also work in 2D (dofs `velx`, `vely`, `rotz`): with
`element_beam_direction_z 0 0 -1` the local y axis flips; with a z vector
in the plane (for example `0 1 0` for a beam along x) the in-plane
bending uses `Iyy`.

Without any of these records the beam is the classic GNU beam: in 2D the
local z axis is the global z axis; in 3D the beam bends only in the plane
of `group_beam_plane` (default the x-y plane) with `Izz`.

## Example

`group_beam_direction_z1.dat` (validation-suite/test-2026): two 3D
cantilevers of length 1 along x, `group_beam_inertia 0 1 4 1`
(Iyy 1, Izz 4), tip force 0.01 in y:

```
group_beam_direction_z 0 0 0 1
element_beam_direction_z 2 0 1 0
```

Element 1 bends about its local z axis (global z) with Izz: tip vely
8.333e-4. Element 2 has local z = global y, so the same force bends it
about its local y axis with Iyy: tip vely 3.333e-3,
`element_beam_force_moment 2 0 0 -0.01 0 0.01 0 0 0 0.01 0 0 0`.

`group_beam_direction_z_reference_point1.dat`: the same beams with
reference points (0.5, 3, 0) (local z = global y) and (0.5, 2, -3)
(local z = -global z).

## Differences with the Professional version

- 3D beam without a local z record: the GNU keeps its planar beam in the
  `group_beam_plane` plane; the Professional uses a full 3D beam with
  its default local z (global z projected perpendicular to the beam,
  global y for a beam along z). Same results for beams in the x-y plane
  under in-plane loads. Owner decision pending (C37,
  [DIVERGENCES-PROFESSIONAL](../DIVERGENCES-PROFESSIONAL.md)).
- Default `group_beam_memory`: GNU `-updated` (local axes follow the
  deformed beam), Professional `-total_linear`; give
  `group_beam_memory -total_linear` to get the Professional results
  (BEAM-MEMORY-DEFAULT).
- `element_beam_direction` has 2 values in the GNU (in-plane direction of
  the planar beam) and is not written for a beam with local axes; the
  Professional writes 9 values (BEAM-DIRECTION-LENGTH).
- A beam whose axial displacement is free stops the GNU direct solver;
  the Professional solves it (BEAM-AXIAL-FREE).
- The perpendicularity check of `group_beam_direction_z` /
  `element_beam_direction_z` uses a tolerance of 1e-3 on the cosine
  (the Professional tolerance is not known).
