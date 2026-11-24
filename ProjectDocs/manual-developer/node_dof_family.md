# node_dof family — spurious `node_dof` requests (refinement/split + eigen)

## Summary

Two independent leftovers/regressions requested the `node_dof` record on
paths where that record legitimately does not exist, aborting with
`Error detected for data item : node_dof`. Both are now fixed:

1. `top.cc` `top()` — a discarded debug call at input-parse time.
2. `so.cc` `solve()` — the band-solver fill-dofs branch running during the
   `control_eigen` task after `VERSION_NEW` had been deleted.

Corpus sfnet 197: **167 PASS / 30 RUNFAIL / 0 PARSE -> 173 PASS / 24 RUNFAIL
/ 0 PARSE** (6 RUNFAIL -> PASS, 0 regressions). Internal suite 16/16.

## Call site 1 — `top.cc` debug leftover (refinement/split family)

### Symptom

`examp7`, `refine4`, `examp15` (pure mesh models: no material, no condif,
therefore `nuknwn == 0` and no `node_dof` record at all) abort during input
processing with:

```
Error detected for data item : node_dof
Error detected for record : 1
```

Backtrace (instrumented `db_error`, `addr2line`): `db_error` <- `db_dbl` <-
`top()`.

### Root cause

`top()` carried a debug block (introduced in `6ab35cb`, the variable removed
in `db6876c` but the call left behind):

```cpp
// debug: dump node 1 stress dofs after the input parse
if ( db_active_index( NODE, 1, VERSION_NORMAL ) ) {
  db_dbl( NODE_DOF, 1, VERSION_NORMAL );
}
```

It runs for **every** model that has node 1, unconditionally requesting
`NODE_DOF` record 1. Any model whose node 1 has no dof (mesh-only models, and
other models where the dof is created later / is inactive at parse time)
aborts. The returned pointer was never used, so the block was dead code with
the only side effect of marking `NODE_DOF` as read.

### Fix

Delete the block (`top.cc`, in `top()`). `NODE_DOF` is read throughout the
solve anyway, so `db_read[]` bookkeeping is unaffected.

## Call site 2 — `so.cc` eigen fill-dofs (eigen family)

### Symptom

`eigen1`, `eigen3`, `ho_heat1` (models WITH dofs) abort after the band solve
prints `number of equations/size of matrix= ...`, with the same
`node_dof` error at record 1. Backtrace: `db_error` <- `db_dbl` <- `solve()`
<- `step_close()` <- `top()`.

### Root cause

`control_eigen` is invoked from `step_close()` (`top.cc:1201`), which runs
**after** the end-of-step sequence in `top()`:

```cpp
db_version_copy( VERSION_NEW, VERSION_NORMAL );   // top.cc:612
db_version_delete( VERSION_NEW );                 // top.cc:613
step_close( YES, ... );                           // top.cc:619 -> solve(-CONTROL_EIGEN)
```

So during `solve(-CONTROL_EIGEN)`, `NODE_DOF` at `VERSION_NEW` does not
exist. The eigen task sets only `band_solver = 1` (`so.cc:462`). Commit
`5b7bb8b` ("fix band fill dofs") added `|| band_solver` to the fill-dofs
condition in `solve()`:

```cpp
if ( bicg_solver || petsc_solver || superlu_solver || lapack_solver
     || band_solver )
  node_dof_new = db_dbl( NODE_DOF, inod, VERSION_NEW );
```

That made the eigen task enter a branch that reads and writes
`NODE_DOF[VERSION_NEW]` -> `db_error`.

The original Tochnog jan-2014 source
(`external-downloads/sfnet/extracted-sources/Sources/Tochnog-Latest-jan-2014/so.cc:715`)
deliberately did **not** include `band_solver` in that condition:

```cpp
if ( bicg_solver || petsc_solver || superlu_solver || lapack_solver)
  node_dof_new = db_dbl( NODE_DOF, inod, VERSION_NEW );
```

`-matrix_lapack` already sets `lapack_solver = 1` (plus `band_solver = 1`),
so the `|| band_solver` addition was redundant for the band/LAPACK route and
only affected the eigen task, where `band_solver` is the only flag set.

Note: `control_eigen` / `node_eigen` do **not** exist in the Tochnog
Professional binaries (neither 25-10-2023 nor 02-08-2026 contain the strings
`control_eigen` / `node_eigen`), so the Professional is not a reference for
this path; the reference is the original GNU semantics.

### Fix

Exclude the eigen task from both the read and the write of the fill-dofs
block in `so.cc` `solve()`. The write block must use the same guard, otherwise
`node_dof_new` would be left NULL:

```cpp
if ( task!=-CONTROL_EIGEN &&
     ( bicg_solver || petsc_solver || superlu_solver || lapack_solver
       || band_solver ) )
  node_dof_new = db_dbl( NODE_DOF, inod, VERSION_NEW );
...
if ( task!=-CONTROL_EIGEN &&
     ( bicg_solver || petsc_solver || superlu_solver
       || lapack_solver || band_solver ) ) {
  tmp = solve_b[ilocal];
  ...
  node_dof_new[iuknwn] += tmp;
}
```

This restores the original eigen behaviour (only `NODE_EIGEN` is filled) while
keeping the band fill-dofs for every non-eigen solve.

## Files and lines changed

- `top.cc` — `top()`: removed the 4-line debug block after the `TIME_CURRENT`
  initialisation (was lines ~244-247).
- `so.cc` — `solve()`: `task!=-CONTROL_EIGEN &&` added to the two fill-dofs
  conditions (was lines ~1147 and ~1158).

## Verification

- Reproduced all 6 failures on the 84988af baseline binary (`Error detected
  for data item : node_dof`).
- After the fix all 6 run rc=0 and reach their targets:
  - `examp7`: 256 elements / 289 nodes = Professional 25-10-2023.
  - `refine4`: 64 elements / 41 nodes = Professional.
  - `examp15`: node 545 = (0.96875, 0.96875) = Professional; 913 elements /
    507 nodes = Professional.
  - `eigen1` / `eigen3`: `control_eigen_values[0] = 0.015408874649`
    (analytic 1D FEM value 0.01540887).
  - `ho_heat1`: `control_eigen_values[0] = 0.00763666044548`, matching the
    exact bilinear-quad 16x16 FEM stiffness eigenvalue to all digits
    (independent numpy check). The test target `0.00768589 +/- 1e-4` is a
    slightly stale Lobatto-era value; the computed value is the exact FEM
    value and passes the tolerance.
- Corpus sfnet 197 A/B (clean builds from 84988af, same runner, LC_ALL=C):
  167/30/0 -> 173/24/0; the only class changes are the 6 target tests
  RUNFAIL -> PASS, 0 regressions.
- Internal suite `scripts/build_safe.sh --clean`: 16/16.

## Evidence

`/tmp/opencode/nodedof-evidence-20260914-101024/` — repro runs, backtraces,
Pro comparisons, corpus classifications, per-test validation.

## Pending / not in scope

- `control_eigen` remains GNU-only (no Professional equivalent); the fix
  restores the original GNU semantics, not a Professional alignment.
- `ho_heat1`'s target is stale (0.64 % high) but within its own 1e-4
  tolerance; not recalibrated here to keep the change narrow.
