# element_beam_force_moment (developer)

## Where

- `database.cc`: `ELEMENT_BEAM_FORCE_MOMENT`, DOUBLE, length 12,
  `version_all`, class ELEMENT (next to `ELEMENT_BEAM_MOMENT`). Not
  print_only: it is the state of `beam_local()` and must survive a
  restart.
- `top.cc`: `db_allocate( ELEMENT_BEAM_FORCE_MOMENT, max_element,
  VERSION_NEW, MINIMAL )` under `any_beam`, next to `ELEMENT_BEAM_MOMENT`
  (no allocation inside the parallel element loop).
- `beam.cc`, `beam_local()`: local end forces `f = f_old + K T du`
  (see [group_beam_direction_z](group_beam_direction_z.md)); PUT
  VERSION_NEW.
- `beam.cc`, `beam_2d()` (classic planar beam, also the 3D
  `group_beam_plane` path): after `ELEMENT_BEAM_MOMENT` is PUT, the plane
  forces are rotated with the rows of `new_rotation_matrix`:
  force_x = row 0 . (f0, f1), force_y = row 1 . (f0, f1), moment_z = m,
  per node; the other 6 values are 0. Output only (the classic beam keeps
  `ELEMENT_BEAM_MOMENT` as its state).
- `check.cc`: needs materi_displacement or materi_velocity_integrated,
  beam_rotation, ndim 2..3.
- `print_beam_force_moment.cc`: used when `ELEMENT_BEAM_MOMENT` is
  absent (local-axes beams).

## Implementation notes

- Sign: the values are the element end forces `K u` (forces of the nodes
  on the beam), which gives the manual minus at the first node (manual
  example `0 0 0 0 0 -10 0 0 0 0 0 10` = Professional probe `cmom`).
- force_x is zero up to round-off (no axial stiffness); no snapping.

## Verification

`element_beam_force_moment1` (manual example on the classic beam),
`element_beam_force_moment2` and the `group_beam_direction_z*` /
`element_beam_direction_z1` tests = Professional 02-08-2026 within 1e-5;
`element_beam_force_moment_2d1` analytic (the Professional 2D beam is
singular on it). The 15 existing beam inputs (test-2014 `bmom*`,
`genbeam1`, test-2026, sfnet `trubea1-4`, `trussbea`, `ho_othr1`,
`tslv_bsym`) give the same `.dbs` and print files as before apart from the
new record.

## Pending refactorings / open items

- `control_print_element -element_beam_force_moment` (manual 6.285) is
  not implemented (`print_el.cc` reads `ELEMENT_BEAM_MOMENT`).
- BEAM-AXIAL-FREE, BEAM-MEMORY-DEFAULT, C37 (see the divergence index).

## Differences with the Professional version

See the user page: BEAM-AXIAL-FREE, BEAM-MEMORY-DEFAULT, C37; 2D not
measured.
