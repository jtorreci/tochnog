# materi_plasti_visco_apply — developer

## Implementation

- **Enum**: `MATERI_PLASTI_VISCO_APPLY` (`tochnog.h` / `tochnog-mod.h`),
  after `MATERI_PLASTI_TENSION_APPLY`.
- **Record** (`database.cc`): INTEGER, length 1, `no_index 1`, CONTROL.
- **Gate**: `set_stress()` (`stress.cc`) sets `viscoplasti = 0` when
  `materi_apply_gate_off( MATERI_PLASTI_VISCO_APPLY,
  CONTROL_MATERI_PLASTI_VISCO_APPLY )` (replaces the control-only
  `control_materi_gate_off`; same result without a global record).

## Verification (Professional binary 02-08-2026 as oracle)

1D bar E 2, tension 1, visco_power 0.1 1, two steps (C33 test model):
with visco 2.666667 / eppxx 0.166667; global -no 1 / 1; control 1 -no 1 / 1;
global -no + control 1 -yes 2.666667 / 0.166667; `group_plasti_apply 0 -no`
3 / 0. GNU equal in every case.
Tests: `materi_plasti_visco_apply1`, `control_materi_plasti_visco_apply1`.

## Pending

None.
