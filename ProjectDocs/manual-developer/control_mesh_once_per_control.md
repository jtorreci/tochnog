# One-shot mesh operations once per control index — developer notes (F4 T7)

## Where it is implemented

- `top.cc`: file-static `step_start_control_entry`, set in `top()` right
  after the `ICONTROL` record is written for a control index (also on a
  `control_repeat` re-entry through `repeat_point`).
- `top.cc` `step_start()`: reads and clears the flag into
  `first_of_control`; `macro()`, `mesh_switch()`, `mesh_mirror()` and
  `mesh_copy()` run only when it is set.

## Implementation details

`step_start(YES, ...)` is called on EVERY time step of a
`control_timestep` index (and `step_start(NO, ...)` once for an index
without time stepping). The mesh operations at the top of `step_start()`
read their control record at the current `icontrol`, so with the record
and `control_timestep` on the same index they ran on every step.
`macro()` has no guard of its own (unlike the interface generation in
`generate.cc`, which skips elements whose `ELEMENT_MACRO_GENERATE` equals
the current index).

A step that is retried after a failed iteration (`goto start_of_step`)
does not see the flag again: the mesh change lives in `VERSION_NORMAL`
and is not rolled back, so it must not be re-applied.

Professional probes (02-08-2026, omen, F4 T7): a 3x2-node mesh with the
operation and a two- or three-step `control_timestep` on index 1:
macro 12 nodes (GNU before 24), copy 12 (GNU 48), mirror 12 (GNU 48),
switch applied once with two steps (GNU twice = identity).
`control_mesh_move` had no effect in the Professional (with or without
time stepping), so it is NOT guarded; `control_mesh_rotate_angle` is only
a parameter of `control_mesh_rotate` in the Professional, and the 3D
`control_mesh_rotate` probe did not run on either binary — both left
unguarded (BACKLOG C46).

Tests: `mesh_macro_timestep1` (202 nodes, = Pro; old 4040) and
`mesh_switch_timestep1` (node 3 `0 0.5`; old `0.5 0`), with node checks
on the `.dbs` in `scripts/build_safe.sh`.

## External features / libraries

None.

## Hardcoded parameters

None.

## Pending refactorings

- The other `step_start()` mesh operations (`adjust_geom`, move, rotate,
  delete/keep, remove) still run on every step; they are idempotent or
  unconfirmed against the Professional (BACKLOG C46).
- A per-operation registry (operation, once/always) would replace the
  scattered `first_of_control &&` tests.
