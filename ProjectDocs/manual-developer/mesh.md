# mesh (short alias of options_mesh)

## Implementation

- **File**: database.cc, `db_number()`.
- The exact-name loop would return the dead `MESH` placeholder (no
  type/class registered, so `no_index=0`): the parser would read
  `-fixed_in_space` as an illegal index. The keyword is intercepted at
  the top of `db_number()` and routed to `OPTIONS_MESH`:
  ```c
  if ( !strcmp( str, "mesh" ) )
    return OPTIONS_MESH;
  ```
  `OPTIONS_MESH` is INTEGER, length ndim, `no_index=1` - exactly the
  Professional record. (The same routing exists in the second alias
  chain of db_number; the first match wins.)

## Verification

- slope_classical_numerical parses `mesh -fixed_in_space
  -fixed_in_space` and runs (the collapse-time target is not met:
  physics, see SEGUIMIENTO).

## Differences with the Professional version

- **default `mesh` kinematics** (known difference, not decided; BACKLOG C8).
  GNU: `-follow_material` (nodes move with `materi_velocity` and no
  `materi_displacement`). Professional 02-08-2026: `-fixed_in_space` (manual
  6.806); 1D bar +50% `-updated`: Pro 0.5, GNU 0.40714. Status: Open,
  MEDIUM, effort large (corpus relies on the moving mesh); tied to C18/C19.
  Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
