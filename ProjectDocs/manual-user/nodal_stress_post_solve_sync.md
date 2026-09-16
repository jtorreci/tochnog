# nodal_stress_post_solve_sync

Status: **active since 2026-09-15**.

## What it does

When a model uses the nodal stress field (`materi_stress`), the stress
reported at a node is re-evaluated from the **solved** velocity at the end of
the time step (after the equilibrium iterations have finished). The reported
stress is therefore in sync with the final solved velocity instead of trailing
it by one iteration.

Only the final state of a step is reported, so the re-evaluation runs once per
step (since 2026-09-16); earlier it ran once per equilibrium iteration, which
cost an extra element pass on every iteration for no observable benefit.

This matches the Tochnog Professional convention: the Professional reports the
nodal stress consistently with the converged (or current) velocity in the same
pass.

## Why it is useful

The stress dofs are not part of the solved unknown vector; historically they
were updated after the solve from the constitutive stress computed with the
**pre-solve** velocity. As a result the reported stress lagged the velocity by
one iteration. The effect is only visible when the equilibrium loop runs a
small number of iterations (a partially deleted element with a reset stress
and a non-conforming tie is the canonical case: `mpc5.dat`).

After this change, `mpc5.dat` reports `sigyy = 0.5` at one iteration, exactly
like the Professional, instead of `1.0`.

## Input syntax

There is no keyword: the behaviour is automatic whenever the absolute nodal
stress recovery is active (i.e. `materi_stress` without a convected or
independent total-strain formulation). No input change is required.

## Parameters

None. The pass covers the six stress components in the standard compacted
order (`sigxx, sigxy, sigxz, sigyy, sigyz, sigzz`).

## Minimal example

`mpc5.dat` (the regression test) is the minimal case:

```
materi_velocity
materi_displacement
materi_stress
...
control_reset_dof   15  -sigyy
control_reset_value_constant  15  1.
control_timestep    20  1. 1.
control_mesh_delete_geometry 20  -world
control_mesh_delete_geometry_element 20  -quad9
control_mesh_delete_geometry_factor 20  0. 0.5
mpc_linear_quadratic -yes
post_point 10  1. 2.5
target_item  0  -post_point_dof 10 -sigyy
target_value 0  0.5 1.e-3
```

With `control_timestep_iterations 20 1` the run now reports `sigyy = 0.5`
(was `1.0`) and passes the target.

## Scope

The re-evaluation applies only to the absolute weighted-average recovery. A
convected stress field (`materi_strain_total`) keeps the historic incremental
update and is not re-evaluated.
