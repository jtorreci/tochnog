# timestep_predict_velocity / timestep_iterations_automatic_apply — developer notes

Sprint 12 lot 3 (manual Professional 6.1089-6.1090). Enums
`TIMESTEP_PREDICT_VELOCITY`, `TIMESTEP_ITERATIONS_AUTOMATIC_APPLY`
(plain records, `no_index = 1`).

## timestep_predict_velocity - PARTIAL (documented)

The Professional's semantics: the initial guess of the linear solve
becomes the previous timestep velocities (instead of zero) when there
is no inertia and no convection. The GNU's staggered scheme CANNOT map
this: `so_bicg.cc` always starts the iteration from x = 0 (r1 = b*p)
and `so.cc` ADDS the solution to the unknowns
(`node_dof_new[iuknwn] += tmp`). Starting from the previous
velocities would double-count them. The record is registered and
behavior-neutral (verified `tslv_part`).

## timestep_iterations_automatic_apply - wired

In `top.cc`, declared per control block, read once before the
automatic branch:

```c
timestep_iterations_automatic_apply = -YES;
db( TIMESTEP_ITERATIONS_AUTOMATIC_APPLY, 0, ..., GET_IF_EXISTS );
if ( timestep_iterations_automatic_apply!=-NO &&
     db_active_index( CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC, icontrol, ... ) )
```

With `-no` the `control_timestep_iterations_automatic` record is not
read and the plain `control_timestep_iterations` / default-2 branches
apply (`use_control_timestep_iterations_automatic` stays 0).

## Gotchas

- The GNU's quasistatic staggered collapses to its fixed point in one
  outer pass, so iteration-count differences are INVISIBLE in these
  models (measured: `control_timestep_iterations 12` gives bit-identical
  results as 2). The automatic machinery only shows in
  residue/dynamic configurations - that is why the gate is verified by
  construction plus the unchanged-physics test (`tslv_auto`) rather
  than by an A/B on iteration counts.
