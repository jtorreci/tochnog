# force_gravity_geometry — developer

## Files and functions

- `database.cc` — `FORCE_GRAVITY_GEOMETRY` (INTEGER, 2, `no_index`, class
  FORCE_GRAVITY, requires `FORCE_GRAVITY`); enum after `FORCE_GRAVITY_TIME`
  in `tochnog.h` / `tochnog-mod.h`.
- `force.cc` — `force_gravity_element_inside(element)`: 1 without the
  record; otherwise 1 only when `geometry()` (`NODE_START_REFINED`,
  `PROJECT_EXACT`) finds every node of the `ELEMENT` record inside.
- `materi.cc` — `materi()`: after the `mesh_activate_gravity` factor, the
  gravity vector is zeroed for an element outside.
- `truss.cc` — `truss()`: the same after `force_gravity_calculate`.

## Implementation details

`force_gravity_calculate()` itself is unchanged: it is also used for nodal
masses (`dof.cc`), water pressures (`groundfl.cc`, `groundda.cc`,
`bounda.cc`, `area.cc`) and `post_element_force`, which the manual does not
restrict.

Measured on the Professional 02-08-2026 (two quad4 columns 0..1 and 1..2,
a quadrilateral -0.5..x): x = 0.9 none, x = 1.0001 .. 1.9 the first column
only (also from 1.55 on, with the centre of the second inside), x = 2.5
both; two vertical trusses, the one outside unloaded. A first sweep with
the quadrilateral points in ring order (a crossed quadrilateral) gave
erratic results: the order is the quad4 one (manual 6.544).

## History

- F5 T5 U5 (2026-11-18): added.

## Hardcoded parameters / pending refactorings

- The check runs per integration point in `materi()`; a per-element cache
  would avoid repeating `geometry()` on large meshes.
- Beam and spring elements are not restricted (not measured).
