# group_spring_direction — developer

Record GROUP_SPRING_DIRECTION (database.cc, DOUBLE, ndim values), consumed
in `spring()` (`spring.cc`).

## Kinematics (C88, F4 T6f, 2026-10-19)

- Direction: the record, normalized; without it (`-spring2` only) the
  INITIAL axis `initial_coord[1] - initial_coord[0]` (before C88: the
  current axis `new_coord`). `initial_coord` = node coordinates, minus
  the integrated displacement `veli` for a `-follow_material` mesh.
- Elongations are projections on that direction: `-spring2`
  `old/new/initial_length = dir . (x2 - x1)` of the old / new / initial
  coordinates, `-spring1` `dir . (x - x_initial)`. `incremental_length =
  new_length - old_length`; the nonlinear diagram uses `new_length -
  initial_length` as total strain (unchanged formulas, new lengths).
  Before C88 the lengths were `array_distance` values (full motion,
  unsigned for spring1).
- Nodal force (element_rhside): `-spring1` `-F dir`, `-spring2` `+F dir`
  at node 0 and `-F dir` at node 1. Tangent: `k dtime dir_i dir_j`
  (signed), `+` on the same node, `-` between the two nodes of a
  spring2; element_lhside gets the diagonal. Before C88 the tangent
  used `|dir_i dir_j|`, and the `-spring1` branch wrote the rhside and
  the matrix at stale `indx` / `jpuknwn` with a negative stiffness.
- `ELEMENT_SPRING_DIRECTION` stores the direction used, `ELEMENT_SPRING_FORCE`
  the force (VERSION_NEW).
- `group_spring_memory` is not consumed: the Professional gives the same
  projection for every memory value (rotation probe force 0).

## Verification

Professional 02-08-2026 probes (2D, k 10): oblique spring2 along (1,1)
moved (0.3,0.1) 2.828427124746 in 1 and 10 steps; no-record spring2
rotated (0,1) 0 for the default, `-total_linear`, `-updated_linear` and
`-updated`; spring1 along (1,0) moved (-0.2,0.5) -2 with node_rhside +2;
two oblique springs on a free node = analytic. Tests (test-2026):
`spring2_oblique_projection1`, `spring2_rotation_initial_axis1`,
`spring1_projection1`, `spring_oblique_free1`, `spring2_updated_rotation1`,
and the C86 section tests `post_element_force_2d_spring1_interior1`,
`post_element_force_3d_truss_spring1`.

## Pending

- The contact springs (`conspr.cc`) are a separate element, not checked.
