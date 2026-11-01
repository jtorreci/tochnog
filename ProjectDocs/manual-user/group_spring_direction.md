# group_spring_direction

Direction of a spring (manual Professional 6.764): `group_spring_direction
index dir_x dir_y [dir_z]`, the index is the element group. Required for a
`-spring1`; for a `-spring2` without the record the direction is the axis
from the first to the second node in the INITIAL geometry.

The spring works only along this direction:

- elongation of a `-spring2` = `dir . (u2 - u1)`, of a `-spring1` =
  `dir . u` (u = displacement of the node from its initial position);
- spring force = `group_spring_stiffness` x elongation (or the
  `group_spring_stiffness_nonlinear` diagram), limited by
  `group_spring_plasti`;
- the force acts along the direction: `+F dir` on the first node and
  `-F dir` on the second node of a `-spring2`, `-F dir` on the node of a
  `-spring1` (a positive force is tension).

Motion perpendicular to the direction gives no force, and the direction
does not rotate with the nodes, so the spring is linear: one step or many
give the same force.

Example (an oblique spring):
```
element 1 -spring2 1 2
group_type 0 -spring
group_spring_stiffness 0 10.
group_spring_direction 0 1. 1.
```
Moving node 2 by (0.3, 0.1) gives the force 10 x 0.4 / sqrt(2) = 2.828.

## Differences with the Professional version

None measured (F4 T6f C88, 2026-10-19: oblique and rotated `-spring2`,
`-spring1`, two oblique springs on a free node, every memory value).
Before C88 the GNU took the elongation as the change of the node
distance, with a tangent `k |dir_i dir_j|`: wrong forces for any spring
moving off its direction, step-dependent results, and a singular matrix
for a (1,-1) spring.
