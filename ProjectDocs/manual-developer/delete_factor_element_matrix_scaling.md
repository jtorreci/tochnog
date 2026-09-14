# delete_factor_element_matrix_scaling

## Files and functions

- `elem.cc` — `void elem( long int element, long int ithread )`.
  - Deleted-element branch (around line 352): the old diagonal-only restore
    was removed; only the nonlocal-softvar early return is kept.
  - `add_to_system_vectors:` label (around line 943): the full-matrix scaling
    and the frozen released-stress RHS substitution are applied here.

## Implementation details

When an element is partially or fully removed, the Professional stores:

- `element_delete_factor` — the remaining stiffness/force fraction.
- `element_rhside_delete` — the frozen released-stress internal force.

The Professional does **not** store an `element_matrix_delete` record. The
deleted element keeps contributing its **full** recomputed element matrix
scaled by `element_delete_factor`, with the frozen released-stress RHS.

The old GNU code (`elem.cc:352`) instead restored a **diagonal-only**
snapshot of the element matrix and the frozen RHS, then jumped directly to
`add_to_system_vectors`. Because `element_matrix` is zeroed at the start of
`elem()`, a deleted element contributed only a lumped diagonal spring plus
the frozen RHS. That under-relaxed the released load and slowed the
staggered u-sigma fixed point (delete3 needed ~6 iterations where the
Professional needs 1).

The new flow:

1. The early diagonal-snapshot branch is gone; a deleted element runs the
   normal element computation and assembles its full consistent matrix.
2. At `add_to_system_vectors`, if `ELEMENT_DELETE_FACTOR` is active:
   - on the first iteration in which the factor is active, the freshly
     computed `element_rhside` is stored as `ELEMENT_RHSIDE_DELETE` (the
     frozen released-stress RHS);
   - on later iterations the frozen RHS is reused;
   - `element_rhside[i] = element_delete_factor * element_rhside_delete[i]`;
   - the **full** `element_matrix` (`nnol*npuknwn` squared) is scaled by
     `element_delete_factor`.

This matches the Professional `.dbs` layout exactly (factor + frozen RHS, no
matrix-delete record) and is applied to every delete producer:
`control_mesh_delete_geometry_factor` (`delete.cc`) and the failure-based
delete (`failure.cc`, `group_materi_failure_crunching` /
`group_materi_failure_void_fraction`).

## External dependencies

None. Uses the internal database API (`db`, `db_active_index`).

## Hardcoded parameters / pending refactorings

- `ELEMENT_MATRIX_DELETE` and its allocation in `top.cc` /
  `database.cc` remain registered but are now unused (the Professional has
  no equivalent record). They can be retired in a later cleanup, together
  with the enum slot, if a renumbering of the database enum is acceptable.
- The frozen RHS is stored for the whole analysis the first time the factor
  becomes active. This mirrors the Professional; the general time-dependent
  softening case (factor history) is not exercised by the corpus.
- The mpc5 staggered lag (released stress routed through the non-conforming
  tie) is a separate mechanism and is NOT closed by this scaling; see the
  SEGUIMIENTO row and the convergence notes.
