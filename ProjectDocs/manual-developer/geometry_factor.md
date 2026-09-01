# geometry_factor — developer notes

Professional manual 6.527. Alias only — the record itself is the GNU
`GEOMETRY_BOUNDA_FACTOR` (DOUBLE, `data_length = 3`, variable length,
`data_class = GEOMETRY`).

## Implementation

One line in the Professional-name translation chain of
`db_number()` (`database.cc`):

```c
else if ( !strcmp( str, "geometry_factor" ) )
  return GEOMETRY_BOUNDA_FACTOR;
```

The chain is reached only after the exact-name loop, so the GNU name
`geometry_bounda_factor` keeps working. Both the keyword detection and
the end-of-variable-values detection in `input.cc` (which checks
`db_number(str) >= 0`) see the translated item.

## How the record is consumed

`geometry.cc::geometry()` evaluates the factor when a node projects on
the entity with the matching index:

- line, `length == 2`: linear — `f[0]*(1-xi) + f[1]*xi` (xi = local
  line coordinate of the projection);
- line, `length == 3`: parabolic — Lagrange quadratic in `xi` mapped to
  [-1,1] with values at start/middle/end;
- triangle, `length == 3`: linear in the three corner shape functions
  (`l0,l1,l2`).

The boundary/force routines that call `geometry()` therefore get the
scaled load for free (bounda.cc, force.cc).

## Verification

- `matrix2` (condif_temperature): linear factors 1..4 per side, middle
  point target temp = 2.5. GNU rc=0 (`post_point_dof` = 2.4999998123)
  against the Professional binary 25-10-2023 (`2.500000000000e+00`).
- `temp2`, `matrix4`: same pattern, rc=0. Corpus 118 -> 121 PASS.

## Gotchas

- The manual's example (line 1..4, node at x=0.2 -> 20*1.6) matches the
  GNU linear interpolation exactly — no semantics difference, only the
  name.
- `control_mesh_delete_geometry_factor` (used by `delete2`/`delete3`) is
  a separate record of class CONTROL; it is NOT implemented yet
  (`delete3` stays RUNFAIL: sigyy 0 vs target 0.6, factor not applied).
