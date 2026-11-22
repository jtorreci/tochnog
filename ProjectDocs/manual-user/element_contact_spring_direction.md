# element_contact_spring_direction

`element_contact_spring_direction index dirNx dirNy dirNz dirT1x dirT1y dirT1z dirT2x dirT2y dirT2z`
(manual 6.418)

After the calculation, this record holds the directions that each contact spring element used. The
GNU name of the record is `element_contactspring_direction`; see
[names_2026_aliases](names_2026_aliases.md).

How the directions are built (the same as the Professional 02-08-2026):

- **Normal N.**
  - With `group_contact_spring_direction`, the given vector is normalised, so `0 2 0` acts as
    `0 1 0`.
  - The normal is then turned, if needed, to point from the side of the first node to the side of
    the second one. A side is the centroid of an element of that node that is not a contact spring.
    This makes FN positive when the spring extends (manual 6.603), whichever sign you give. For
    example, `0 0 1` and `0 0 -1` give the same normal.
  - Without a direction, or with `group_contact_spring_direction_automatic -yes`, the normal is
    that side direction itself; see
    [group_contact_spring_direction_automatic_planes](group_contact_spring_direction_automatic_planes.md).
- **First tangent T1** = e_z x N, normalised: in 2D, T1 = (-N_y, N_x). When N lies along z,
  T1 = e_y x N instead, so N (0,0,1) gives T1 (1,0,0).
- **Second tangent T2** = N x T1.

Example: conspr6 with `group_contact_spring_direction 1 1 1 1` gives N (-0.57735, -0.57735, -0.57735),
T1 (0.707107, -0.707107, 0) and T2 (-0.408248, -0.408248, 0.816497).

Tests: `contact_spring_3d_direction1`, `contact_spring_3d_direction2`, `contact_spring_2d_direction1`.

## Differences with the Professional version

When the normal force FN is exactly 0, the GNU version follows manual 6.608: the tangential force cannot
exceed f FN, so the spring slides with no tangential force. The Professional 02-08-2026 keeps the elastic
tangential force there, even with f = 0.5. That is not physical, so the GNU version does not follow it.
With any nonzero FN, both limit the tangential force to f |FN|. Guard test:
`contact_spring_friction_zero_normal1`.

