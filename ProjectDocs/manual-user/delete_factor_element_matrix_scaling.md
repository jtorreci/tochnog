# delete_factor_element_matrix_scaling

## Description

Controls how a partially deleted element keeps contributing to the
equilibrium when its delete factor is between 0 and 1 (a full removal uses
the plain `delete_element` path instead).

This behavior is reached through:

- `control_mesh_delete_geometry_factor` — progressive removal of the
  elements that fall inside a geometry (see the mesh-delete family).
- `group_materi_failure_crunching` and
  `group_materi_failure_void_fraction` — elements marked for removal by a
  failure criterion.

## Semantics

For an element with `element_delete_factor` `f`:

- the internal force released at the moment the element is first marked
  (the frozen released-stress right-hand side) is scaled by `f` and applied
  as a nodal load from then on;
- the element contributes **nothing else** to the system: its stiffness
  matrix and its lumped nodal-stress weight are zero.

With `f = 1` the element is untouched. With `f = 0` it is fully removed.
Intermediate values let a deleted region relax its load progressively, which
is how the Professional reproduces the `delete3` test (the stress in the
remaining element settles at `f` times the initial stress) and the `mpc5`
test (the released load reaches the tied linear element in the same solve).

## Notes

- The frozen released load is a constant, so the element must not keep a
  stiffness in the solve: a residual stiffness would respond to the strain
  while the internal force is frozen, which breaks the one-pass Newton
  cancellation of the staggered solve and makes the released load relax
  geometrically instead of settling in the test's own iteration count.
- Because the element no longer updates its nodal stress, the stress of a
  node that only belongs to deleted elements stays at its value at the
  moment of deletion; a node shared with a surviving element takes the
  surviving element's stress. This matches the Professional's field.
- For non-conforming interfaces (MPC ties) the released stress of a slave
  node is redistributed to its masters by the tie elimination, in the same
  solve; see the mpc family.
