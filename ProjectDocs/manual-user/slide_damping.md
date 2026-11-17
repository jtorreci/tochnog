# slide_damping

## Description

`slide_damping index damping_n damping_t` (manual 6.1045): viscous normal
and tangential damping of a `slide_geometry` with the same index. Also:

- `control_slide_damping_apply index switch` (manual 6.370): `-no` switches
  the slide_damping records off in that control index; default `-yes`.
- `slide_user index switch` (manual 6.1044): `-yes` calls the user routine
  `user_slide_friction_force` (user.cc) for the maximum friction force.

## Effect on the calculation

- The damping force on a node of the slide is `-damping_n*vn*n -
  damping_t*vt`, with `v` the node velocity, `vn` its normal component and
  `vt` its tangential part. It acts with or without `slide_stiffness`
  (without it the springs are 0), and it is part of `node_slide_force`
  (slot 0 normal, slot 1 tangential size).
- The `slide_plasti_friction` cap applies to the spring force; the damping
  force is added outside the cap.
- `slide_user -yes`: the routine gets the node, the slide index, the normal
  spring force and the tangential displacement size, and returns the
  maximum tangential force (on input the `slide_plasti_friction` cap). The
  distributed routine stops with `Error: routine user_slide_friction_force
  is not programmed.`, as in the Professional; edit user.cc and rebuild.

## Example

```
slide_geometry 10  -bottom
slide_stiffness 10  1.e3 1.e3
slide_damping 10  100. 10.
control_slide_damping_apply 0 -yes
```

Self-checking inputs in `validation-suite/test-2026`: `slide_damping1`
(damping only, slide force = c*v), `slide_damping2`
(`control_slide_damping_apply -no`, Professional values),
`slide_user_error1`.

## Differences with the Professional version

- The Professional 02-08-2026 computes its "damping" force as `c*u/dt`,
  with `u` the TOTAL slide displacement: a spring of stiffness `c/dt` that
  changes with the time step, not a damping. The GNU applies true viscous
  damping `c*v` (owner decision 2026-11-17, DIVERGENCES
  SLIDE-DAMPING-VISCOUS). To reproduce a Professional input, add
  `damping/dt` to `slide_stiffness` instead.
