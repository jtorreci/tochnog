# delete_factor_element_matrix_scaling

## Files and functions

- `elem.cc` — `void elem( long int element, long int ithread )`.
  - Deleted-element branch (around line 352): only the nonlocal-softvar
    early return is kept (the old diagonal-only restore was removed by
    `7687e0a`).
  - `add_to_system_vectors:` label (around line 940): the frozen
    released-stress RHS substitution, the matrix/lhside zeroing and the
    frozen-RHS capture are applied here.

## Implementation details

When an element is partially or fully removed, the Professional stores:

- `element_delete_factor` — the remaining force fraction.
- `element_rhside_delete` — the frozen released-stress internal force.

The Professional does **not** store an `element_matrix_delete` record. The
deleted element applies its frozen released-stress internal force as a
nodal load and contributes nothing else to the system.

The flow (`elem.cc`, `add_to_system_vectors`), when
`ELEMENT_DELETE_FACTOR` is active:

1. On the first iteration in which the factor is active, the freshly
   computed `element_rhside` is stored as `ELEMENT_RHSIDE_DELETE` (the
   frozen released load); on later iterations the frozen record is reused.
2. `element_rhside[i] = element_delete_factor * element_rhside_delete[i]`.
3. `element_lhside[i] = 0` and `element_matrix[i] = 0`.

Steps 2-3 are the energy-consistent (Jacobian-consistent) pairing: the
frozen released load is a **constant**, so its tangent contribution to the
staggered momentum residual is zero. The historic code
(`7687e0a`) instead scaled the full element matrix by the factor
(`element_matrix[i] *= factor`), which left a residual stiffness in the
solve that responded to the strain while the internal force stayed frozen.
For plain elements the staggered map collapses to one pass because
`A = dt*K = dt*B^T*C*B` exactly cancels the constitutive term; the residual
stiffness of the deleted element broke that cancellation, so the released
load relaxed geometrically (mpc5 contracted at rate 0.5, reaching
`sigyy = 0.5` only at ~12 iterations instead of at the test's own 2).

Zeroing `element_lhside` matters too: it is the lumped weight of the nodal
stress projection (`parallel_new_dof_diagonal`, `dof.cc`). With the frozen
stress row RHS (0 at capture) but a non-zero lhside, the deleted element
diluted the surviving elements' stress recovery at a shared node and slowed
the nodal stress even when the velocity was already one-pass.

This matches the Professional `.dbs` layout (`element_delete_factor` +
`element_rhside_delete`, no matrix-delete record) and its measured field.

## Verification (measured)

Against the Professional binary 02-08-2026, `LC_ALL=C`:

- mpc5 (quad9 deleted 50%, tie, reset `-sigyy 1`): `post_point_dof 10
  -sigyy = 0.5` at the test's own 2 iterations, rc=0 (was 0.7, rc=1). The
  full `node_dof` field matches the Professional to ~5e-17.
- delete3 (top quad9 deleted 40%): `post_point_dof 1 -sigyy = 0.6`, rc=0
  (was already 0.6 at 5 iterations; now at 2). The full field matches the
  Professional to ~1e-14.
- The match holds for every delete factor tested (mpc5 with 0.3/0.5/0.7/0.8,
  delete3 with 0.4/0.6/0.8).
- Corpus: sfnet 197 and validation-suite 233 classification byte-identical
  to the 1e31fd6 baseline; internal suite 16/16.

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
- The GNU nodal-stress projection is a post-solve lumped update of the
  pre-solve element residual (`parallel_new_dof_diagonal`), so the reported
  nodal stress trails the solved velocity by one iteration. The Professional
  reports it in sync. This is a general property of the staggered u-sigma
  scheme (not delete-specific): mpc5 reaches the target at the test's own
  2 iterations but not at k=1 (the Professional is one-pass at k=1). Closing
  that one-iteration reporting lag is a separate work unit (post-solve
  constitutive re-evaluation of the nodal stress).
