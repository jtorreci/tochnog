# One-shot mesh operations: once per control index

## What it does

`control_mesh_macro`, `control_mesh_copy`, `control_mesh_mirror` and
`control_mesh_switch` are applied once per control index, on the first
step of that index, also when the same index carries a
`control_timestep` with several steps. Since F4 T7 this matches the
Professional 02-08-2026.

## Why it is useful

Before F4 T7 the GNU re-applied these records on every time step of the
index. A `control_mesh_macro` next to a `control_timestep` added a copy of
the mesh per step (the Professional `mesh_activate1`: 202 nodes became
4040 after 20 steps, and the 1000-step input ran out of memory), a
`control_mesh_copy` doubled the mesh per step and a `control_mesh_switch`
swapped the axes back and forth.

## Input syntax

No new keyword. Put the mesh record and the time stepping on the same
index, as in the Professional examples:

```
control_mesh_macro 10   -rectangle 0 2 101
control_mesh_macro_parameters 10   0.05 0.5 0.1 1.
control_timestep  10  0.05 1.
```

## Parameters

None. The operation runs again when the control loop enters the index
again (a `control_repeat` re-entry).

## Minimal example

`validation-suite/test-2026/mesh_macro_timestep1.dat` (the Professional
`mesh_activate1` with 20 steps: 202 nodes at the end) and
`mesh_switch_timestep1.dat` (two steps, axes swapped once).

## Differences with the Professional version

- `control_mesh_move` is not part of this rule: in the probes the
  Professional did not move the mesh at all, so the once-per-control
  behaviour could not be confirmed; the GNU still moves the mesh on every
  step of a `control_timestep` index (BACKLOG C46, DIVERGENCES C46).
- On one index with both `control_mesh_copy` and `control_mesh_switch` the
  Professional copies first and then switches; the GNU switches first
  (BACKLOG C46).
- `control_mesh_copy` / `control_mesh_mirror` on an index after the first
  time step stop with "Error detected for data item :
  node_dof_start_refined" in the GNU (pre-existing, BACKLOG C45); the
  Professional runs them.
