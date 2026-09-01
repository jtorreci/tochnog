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
