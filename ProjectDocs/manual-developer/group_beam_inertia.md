# group_beam_inertia

## Where

- Registration: `database.cc` — `GROUP_BEAM_INERTIA` variable length
  1..3 (fixed_length 0); the length is not enforced by the parser.
- Consumption: `beam.cc` (`beam()`), which reads up to 3 values and
  uses `I = value[1]` (Izz) when 2+ values are given, `I = value[0]`
  for the legacy single-value form. The beam matrix is the standard
  Euler-Bernoulli axial + 2D bending matrix in (UX, UY, ROTZ).
- Memory: `beam.cc` maps `-total_linear` onto
  `-UPDATED_WITHOUT_ROTATION` (linear kinematics, initial coordinates
  fixed) before the type check.
- 2D rotation alias: `database.cc` `db_number()` — when `ndim==2` the
  names `rotx`/`roty` resolve to the name index whose registered
  basename is `rotz` (the single in-plane rotation unknown created by
  `beam_rotation` in 2D). Required because the Professional inputs
  prescribe `-rotx -roty -rotz` even in 2D (its beam model keeps the
  three rotation unknowns).

## Design notes

The Professional truss_beam is a full 3D frame element; the GNU beam
is 2D. A/B against the Professional binary (probe with
`group_beam_inertia 0 1. 100. 1.`) showed the in-plane (x-y) bending
uses `Izz` — value 2 — so the GNU consumption reads exactly that
slot. trubea2/trubea3 need OUT-OF-plane bending (3D beam) — PENDING
(the GNU beam has no stiffness for the transverse direction, the
system is singular).

## Pending

- Full 3D truss_beam (bending about both transverse axes + torsion) —
  blockers of trubea2/trubea3.
- `group_beam_direction_z*` and shear-deformable options.
