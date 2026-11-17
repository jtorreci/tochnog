# slide_damping — developer

## Files and functions

- `database.cc` — `SLIDE_DAMPING` (DOUBLE, 2, requires `SLIDE_GEOMETRY`),
  `SLIDE_USER` (INTEGER, 1), `CONTROL_SLIDE_DAMPING_APPLY` (INTEGER, 1,
  class CONTROL); enums after `CONTROL_SLIDE_STIFFNESS_APPLY` in
  `tochnog.h` / `tochnog-mod.h`.
- `slide.cc` — `slide()` routes a slide with `slide_stiffness` OR
  `slide_damping` to `slide_spring()`; the damping and the user hook live
  there.
- `user.cc` — `user_slide_friction_force(inod, islide, normal_force,
  tangential_displacement, friction_force)`, distributed as an error stub.

## Implementation details

`slide_spring()`:

1. `control_slide_damping_apply` (control index, default `-YES`) gates the
   `slide_damping` read (`cn`, `ct`). A `-NO` `control_slide_stiffness_apply`
   now zeroes the springs but keeps the damping; without damping it
   returns as before.
2. `kn`, `kt` = 0 when `slide_stiffness` is absent. Without a displacement
   unknown (`materi_displacement` / `materi_velocity_integrated`) a damped
   slide still runs with zero spring displacement.
3. `slide_user -yes`: `user_slide_friction_force` replaces `cap_t` after
   the `slide_plasti_friction` cap.
4. Damping after the spring loop: `vn = v.n`, `vt = v - vn*n`,
   `fd = -ax*(cn*vn*n + ct*vt)` into `node_rhside`, `ax*(cn*|n_i| +
   ct*(1-|n_i|))` into `node_lhside` (velocity unknown: no `dtime` factor,
   unlike the springs).
5. `node_slide_force`: slot 0 `-S*ax - fd.n`, slot 1 minus the size of the
   tangential part of spring + damping force. `node_slide_f` keeps the
   friction-only yield value.

Measured on the Professional 02-08-2026 (probes of `slide_damping1/2`):
force `(k + c/dt)*u` with `u` the integrated (total) displacement, checked
at dt 0.01 and 0.02 and with `slide_damping` alone; `control_slide_damping_apply
-no` gives the undamped `.dbs`; `slide_user -yes` stops with the stub
message. The GNU keeps the viscous law (DIVERGENCES SLIDE-DAMPING-VISCOUS).

## History

- F5 T5 U2 (2026-11-17): new.
