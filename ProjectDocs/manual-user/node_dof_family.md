# node_dof family — models aborting with `Error detected for data item : node_dof`

## What was wrong

Some models aborted with:

```
Error detected for data item : node_dof
Error detected for record : 1
```

even though the model did not (and should not) have a `node_dof` record.
Two situations were affected:

- **Mesh-only models** (no material and no conduction, so no degrees of
  freedom at all), in particular the global h-refinement / split examples
  `examp7`, `refine4`, `examp15`.
- **Eigenvalue analyses** (`control_eigen`) of otherwise valid models:
  `eigen1`, `eigen3`, `ho_heat1`.

## Why it happened

- The mesh-only models were blocked by a leftover debug statement in the
  input-processing routine that always asked for the degrees of freedom of
  node 1. A model with no dofs has no such record, so Tochnog stopped.
- The eigenvalue path ran a bookkeeping step that tries to copy the solution
  back into the node degrees of freedom. That step is valid for a normal
  time-step solve, but the eigenvalue analysis is invoked after the step
  bookkeeping has already been discarded, so the target record no longer
  exists.

In both cases the request was spurious: the code assumed a `node_dof` record
that the model or the analysis step legitimately lacks.

## What changed

- The leftover debug request was removed.
- The eigenvalue analysis no longer tries to write the node degrees of
  freedom (it only produces the eigenvalue / eigenvector results, exactly as
  in the original Tochnog). Normal solves are unaffected.

## What you get

The six affected sfnet tests now run to completion and reach their targets:

| Test | Result |
| --- | --- |
| `examp7` | rc=0, 256 elements / 289 nodes (identical to Tochnog Professional) |
| `refine4` | rc=0, 64 elements / 41 nodes (identical to Tochnog Professional) |
| `examp15` | rc=0, node 545 = 0.96875 (identical to Tochnog Professional) |
| `eigen1` | rc=0, first eigenvalue 0.015408874649 |
| `eigen3` | rc=0, first eigenvalue 0.015408874649 |
| `ho_heat1` | rc=0, first eigenvalue 0.00763666044548 |

Across the whole sfnet corpus the classification moves from
167 PASS / 30 RUNFAIL to 173 PASS / 24 RUNFAIL with no regressions.

## Notes

- `control_eigen` is a GNU-only feature; it does not exist in the Tochnog
  Professional versions, so these two eigenvalue tests have no Professional
  reference. Their results are checked against the analytic finite-element
  eigenvalues instead.
- The `ho_heat1` test target (0.00768589) is slightly higher than the exact
  finite-element eigenvalue (0.00763666); the test passes within its own
  tolerance. The target was left unchanged to keep the fix narrow.
