# materi_plasti_tension_apply — developer

## Implementation

- **Enum**: `MATERI_PLASTI_TENSION_APPLY` in `tochnog.h` / `tochnog-mod.h`
  (in sync), after `CONTROL_MATERI_PLASTI_VISCO_APPLY`.
- **Record** (`database.cc`, `db_initialize`): INTEGER, length 1,
  `no_index 1`, data_class CONTROL.
- **Gate**: `materi_apply_gate_off( MATERI_PLASTI_TENSION_APPLY,
  CONTROL_MATERI_PLASTI_TENSION_APPLY )` (`general.cc`; control record at
  `ICONTROL` overrules the global one, default `-yes`) at three sites:
  - `materi_direct_cutoff()` and the `tension_direct` read of the direct
    laws (`stress.cc`, `ten_present` / `ten_record`; these replaced the older
    control-only `control_materi_gate_off`, same result without a global
    record);
  - the `GROUP_MATERI_PLASTI_TENSION` branch of `plasti_rule()`
    (`plasti.cc`): new; before F5 the tension yield surface was not gated by
    the control form either.

## Verification (Professional binary 02-08-2026 as oracle)

| Case | Professional | GNU |
|------|--------------|-----|
| quad4 tension_direct 1, global -no | sigyy 21.978022, eppyy 0 | same |
| quad4, control 1 -no | 11.989011 / 0.0093 | same |
| quad4, global -no + control 1 -yes | 1 / 0.0193 | same |
| 1D bar tension 1, global -no | sigxx 3, eppxx 0 | same |
| 1D bar, control 1 -no | 3 / 0 | same (GNU before F5: 1 / 1) |
| 1D bar, global -no + control 1 -yes | 1 / 1 | same |
| 1D bar, global -no + control 0 -yes | 3 / 0 | same |

Tests: `materi_plasti_tension_apply1`, `materi_plasti_tension_apply2`,
`control_materi_plasti_tension_apply_yield1` (test-2026, registered in
`scripts/build_safe.sh`).

## Pending

- `control_materi_gate_off()` (control-only helper) is still used by other
  `control_materi_*_apply` gates without a global form.
