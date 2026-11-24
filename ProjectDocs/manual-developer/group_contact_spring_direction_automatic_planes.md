# group_contact_spring_direction_automatic_planes — developer

## Files and functions

- `database.cc` — `GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC_PLANES` (INTEGER,
  MDIM, class CONTACTSPRING, requires GROUP_TYPE), GNU name
  `group_contactspring_direction_automatic_planes`; the Professional name
  `group_contact_spring_direction_automatic_planes` is translated in
  `db_number()`. Enum after `GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC` in
  `tochnog.h` / `tochnog-mod.h`.
- `conspr.cc` — `contactspring()`:
  - automatic direction: after the centroid difference, components with a
    `-NO` switch are zeroed and the vector is normalized; a zero vector
    stops with an error.
  - friction: without a `GROUP_CONTACTSPRING_FRICTION` record (and without
    `_FRICTION_AUTOMATIC`) `friction_unlimited` skips the slip test.
  - 2D: the first tangential direction is `(-n_y, n_x, 0)`; 3D keeps the
    `e_x × n` rule.
  - element matrix: signed `n_i n_j` products (before `scalar_dabs`).

## Implementation details

Measured on the Professional 02-08-2026 with conspr7 (two quad4 blocks, two
contact springs) plus a horizontal load 0.5 on each top node, and a sheared
top block whose centroid direction is (-0.447, -0.894):

| case | Professional | GNU before |
|---|---|---|
| axis direction, sideways load | T = -0.5 / -0.5, disx 24.5 | T = 0 / 0, disx 1.5e12 |
| inclined direction | N 0.1266 / 2.1094, T -0.2533 / 0.2533 | N 0 / 2.236, T 0 |
| planes `-no -yes` | N 0 / 2.0, T -0.4056 / -0.5944 | record unknown |
| planes `-yes -no` | N 0.3696 / 0.6304, T 0 / 2.0 | record unknown |

GNU now matches every value to about 11 digits. The friction default
(f = 0) removed the tangential force; the 2D tangent along z sent the
tangential force to index 2; the absolute value in the matrix gave another
converged split for inclined directions.

## History

- F5 T5 U6 (2026-11-19): planes record; friction default, 2D tangent and
  signed matrix fixed.

## Hardcoded parameters / pending refactorings

- Without `group_contact_spring_direction` and without
  `group_contact_spring_direction_automatic` the Professional stops with
  an error; GNU takes the automatic direction (not changed).
- The 3D tangential directions are not measured.
