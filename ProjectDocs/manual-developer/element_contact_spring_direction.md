# element_contact_spring_direction — developer

The record is an alias of `ELEMENT_CONTACTSPRING_DIRECTION`. `contactspring()` in `conspr.cc` writes
it (`VERSION_NEW`) with 9 values: N, T1 and T2.

## Code

- **`contactspring_side_direction()`** (static) returns centroid(element of node 1) minus
  centroid(element of node 0). For each node it takes the first element in `NODE_ELEMENT` that is
  not a contact spring, and falls back to the spring itself. Both the automatic normal and the
  orientation of a user normal use it.
- **User normal** (`group_contact_spring_direction`, 2D and 3D): it is normalised (a zero vector is
  an input error), then negated when its inner product with the side direction is negative.
- **First tangent:** 2D uses (-n_y, n_x, 0); 3D uses e_z x n. When that cross product vanishes, the
  normal lies along z, and the code uses e_y x n.
- **Second tangent:** n x t1.

## Professional 02-08-2026 measurement (2026-11-23)

The probes were conspr6 and conspr7 with nine normals, the spring node order reversed, and
compression, tension and shear loads. Every direction triad now equals the Professional one.

| Case | Professional | GNU before |
|---|---|---|
| 3D normal 1 1 1 | node 1 velz -0.5 | Kept the vector unnormalised, velz -0.47064 |
| 2D normal 0 2 0 | Force 1 | Force -0.5 |
| User normal sign | Turned to the side direction | Kept the sign given (FN negative in extension for a reversed normal) |
| 3D first tangent | e_z x n | e_x x n |

The tangent rule changes only the components of the tangential force, not the response, because the
tangential stiffness and the friction limit are isotropic.

## Differences with the Professional version: CONSPR-FRICTION-ZERO-NORMAL (deliberate, owner 2026-11-23)

The case is FN = 0 exactly, with a horizontal normal and a vertical load. The GNU follows the manual
(6.608: the total tangential force cannot exceed f FN), so the spring slides with zero tangential
force. The Professional keeps the full elastic tangential force, even with f = 0.5. With any nonzero
FN (down to 4e-10) both limit the force to f |FN|.

The GNU keeps the manual rule; the Professional behaviour is not physical. Guard test
`contact_spring_friction_zero_normal1`, which the Professional fails. See `DIVERGENCES-PROFESSIONAL.md`.
