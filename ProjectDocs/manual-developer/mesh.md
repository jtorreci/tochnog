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

## Default (F3, C8)

- **File**: top.cc, `top()` (block `if ( materi_velocity )` before the
  time loop). The default written into `OPTIONS_MESH` is now
  `-FIXED_IN_SPACE` (manual Professional 6.806, was `-FOLLOW_MATERIAL`);
  an existing record is read over it with `GET_IF_EXISTS`, and
  `materi_displacement` still forces `-FIXED_IN_SPACE`. The record is
  always PUT, so every reader (locate.cc, general.cc, beam.cc, spring.cc,
  truss.cc, conspr.cc) sees the resolved value.
- conspr.cc keeps a local `-FOLLOW_MATERIAL` initialiser before its
  `GET_IF_EXISTS`; it is dead whenever `materi_velocity` is set (the record
  then exists). Pending cleanup: align it with the new default.
- Tests: `validation-suite/test-2026/mesh_default_fixed_in_space1.dat`
  (bar +50%, sigxx 0.5 = Pro; old 0.407136) and
  `mesh_default_fixed_in_space2.dat` (rigid rotation, stresses 0 = Pro;
  old 0.109567).
- Corpus impact: see `odd/tasks/f3-kinematics.md` (tests calibrated on the
  old moving-mesh default; external sfnet targets listed in BACKLOG).

## Verification

- slope_classical_numerical parses `mesh -fixed_in_space
  -fixed_in_space` and runs (the collapse-time target is not met:
  physics, see SEGUIMIENTO).

## Differences with the Professional version

- **default `mesh` kinematics** (C8): RESOLVED in F3 (2026-10-07). The
  default is now `-fixed_in_space`, as the Professional.
- **explicit `-follow_material`: in-step geometry** (known difference,
  documented). GNU evaluates the velocity gradient of a `-follow_material`
  step on the END-of-step geometry (nodes already advanced by the
  iteration), the Professional on the START-of-step geometry. The results
  converge with the number of steps: 2-quad4 column sheared by 1 with
  `-updated_linear`, node 3: one step GNU 0.47222 / Pro 0.5; 100 steps
  GNU 0.52048 / Pro 0.52093. Status: Open, LOW (F3, 2026-10-07). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).

- **Deliberate divergence (decision 2026-10-08):** with an explicit `-follow_material` the GNU builds the step on the end-of-step geometry and the Professional on the start-of-step geometry (2-quad4 column: 0.47222 vs 0.5 in one step, 0.52048 vs 0.52093 in 100 steps). Kept on purpose; see `ProjectDocs/DIVERGENCES-PROFESSIONAL.md` (F3-INSTEP).
