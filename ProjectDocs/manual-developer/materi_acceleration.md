# materi_acceleration

Developer notes for the acceleration dof group.

## Where

- `tochnog.h` / `tochnog-mod.h`: `MATERI_ACCELERATION` dof-type enum
  (inserted before `MATERI_VELOCITY`), `materi_acceleration` flag and
  `acc_indx` index declarations.
- `initia.cc`: global definitions (`materi_acceleration = 0`,
  `acc_indx = -1`).
- `input.cc` (initia part, ~line 500): the `materi_acceleration` keyword
  handler. Allocates `ndim` VECTOR dofs of type `-MATERI_ACCELERATION`,
  starting at `acc_indx`. NOT a principal unknown (like displacement):
  the acceleration rows have no equations.
- `database.cc` `db_initialize()`: basename resolution of the new dof type
  (`accx`, `accy`, `accz` — Professional labels; note the manual spells the
  y/z components `-accy`/`-accz`, mirroring `velx`/`vely`/`velz`).
- `dof.cc` `parallel_new_dof_diagonal()`: per iteration, after the solve,
  `acc_new[idim] = (v_new - v_old)/dtime` (the same block that integrates
  displacements). `v_old` is the step-start velocity (VERSION_NORMAL is not
  touched during the step), so at convergence the record holds the
  acceleration of the step.
- `bounda.cc` (generic assignment branch): a `bounda_dof` whose label
  resolves inside `[acc_indx, acc_indx + ndim*nder)` is intercepted BEFORE
  the generic `node_dof_new[iuknwn] = load` write. It instead binds the
  velocity dof of the same direction: `node_bounded[vel] = 1` and
  `v_new = v_old + a*dt`. The derived record of the bounded node then equals
  the prescribed acceleration automatically (no extra write needed).

## How the labels resolve

`db_initialize()` walks `dof_type` and registers the dynamic basename names
(`accx`...) as `name[idat]` entries with `dof_label[iuknwn] = -idat`. All
INTEGER-record consumers (`bounda_dof`, `target_item`, `print_history`,
`post_point_dof`, `control_reset_dof`, ...) resolve `-accx` through
`db_number()` + `array_member(dof_label,...)`, so no further registration is
needed once the dof group exists.

## Empirical verification (Professional binary 25-10-2023)

- `dynamic3` (1D truss, prescribed `-accx` = 1.234, dt = 0.1): GNU node 1
  displacement 62.317 at t = 10 == Professional exactly (Euler forward
  integration of the velocity bound). rc=0 (target 61.7 ± 1).
- `dynamic4` (node_mass + force): target disx = 1.0 and accx = 2.0 at
  t = 1.0, GNU rc=0.
- node_dof column order follows the initia order (verified: displacement,
  velocity, acceleration -> disx, velx, accx).

## Gotchas

- In 2D the y label is `accy` and in 3D the z label is `accz` (NOT
  `acy`/`acz`) — cross-checked against the dof_label list of the
  Professional manual.
- The momentum equilibrium of a node with a prescribed acceleration is NOT
  assembled: its velocity row is bounded like any Dirichlet dof. The
  inertial reaction of its mass appears in the boundary reaction, matching
  the Professional.
