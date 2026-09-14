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

- the element stiffness matrix is scaled by `f`;
- the internal force released at the moment the element is first marked
  (the frozen released-stress right-hand side) is scaled by `f` and used
  from then on.

With `f = 1` the element is untouched. With `f = 0` it is fully removed.
Intermediate values let a deleted region relax its load progressively, which
is how the Professional reproduces the `delete3` test (the stress in the
remaining element settles at `f` times the initial stress).

## Notes

- The scaling is applied to the full consistent element matrix, not only to
  its diagonal. This is required to match the Professional and to let the
  deleted load relax in the same solve.
- For non-conforming interfaces (MPC ties) the released stress still reaches
  the master through the tie one solve later; that is a separate staggered
  coupling effect and is not addressed here.
