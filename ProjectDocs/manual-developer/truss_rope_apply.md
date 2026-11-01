# truss_rope_apply — developer

## Implementation

- **Enums**: `TRUSS_ROPE_APPLY`, `CONTROL_TRUSS_ROPE_APPLY`
  (`tochnog.h` / `tochnog-mod.h`), after `GROUP_TRUSS_ROPE`.
- **Records** (`database.cc`): INTEGER, length 1, CONTROL;
  `truss_rope_apply` has `no_index 1`.
- **Consumer**: `truss()` (`truss.cc`), the rope clamp
  (`group_truss_rope==-YES && force<0 -> 0`) is skipped when
  `materi_apply_gate_off( TRUSS_ROPE_APPLY, CONTROL_TRUSS_ROPE_APPLY )`.

## Verification

2D truss (A 1, E 100) compressed 0.1 in 10 steps, node 2 rhside velx:
no rope 5.5 (both binaries); rope 0 (both); rope + `truss_rope_apply -no`
GNU 5.5, Professional 0; rope + `control_truss_rope_apply 10 -no` GNU 5.5,
Professional 0; rope + global -no + control 10 -yes 0 (both). The
Professional ignores the switch; the manual (6.1092 / 6.393) decides.
Tests: `truss_rope_apply1`, `control_truss_rope_apply1`, `truss_rope_apply2`.

## Pending

None.
