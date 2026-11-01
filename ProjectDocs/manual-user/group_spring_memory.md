# group_spring_memory

Spring memory model selector (manual Professional 6.765). Value: `-updated_linear`,
`-total_linear` or `-updated`; the index is the element group (see element_group).

The record is accepted and stored. The spring law does not depend on it:
the elongation is always the relative displacement projected on a fixed
direction (see [group_spring_direction](group_spring_direction.md)), so the
spring is linear and gives the same force in one step or in many. This is
what the Professional 02-08-2026 does for every memory value (measured
2026-10-19: a `-spring2` whose axis rotates by 90 degrees carries no force
with `-total_linear`, `-updated_linear` and `-updated` alike). As in the
Professional, `-updated_linear` and `-updated` need
`materi_velocity_integrated`, `-total_linear` needs `materi_displacement`.

Example:
```
group_type 0  -spring
group_spring_memory 0  -total_linear
group_spring_stiffness 0  1.0
```

## Differences with the Professional version

None measured. The manual describes `-updated` as geometrically nonlinear
for large spring rotations; the Professional does not show that behaviour
(the rotated spring above carries no force) and the GNU follows the
measured behaviour.
