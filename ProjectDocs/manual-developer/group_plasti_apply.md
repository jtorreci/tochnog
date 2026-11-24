# group_plasti_apply — developer

## Implementation

- **Enum**: `GROUP_PLASTI_APPLY` (`tochnog.h` / `tochnog-mod.h`), after
  `PLASTI_APPLY`.
- **Record** (`database.cc`): INTEGER, length 1, data_class MATERI,
  `data_required GROUP_TYPE`, indexed by the element group.
- **Consumer**: `set_stress()` (`stress.cc`), right after
  `plasti_off = materi_apply_gate_off( PLASTI_APPLY, CONTROL_PLASTI_APPLY )`:
  a `-no` for the current group `gr` sets `plasti_off = 1`, so every path
  gated by `plasti_apply` (plasti_incr, hypoplasticity, direct cut-offs,
  yield test) is skipped for that group. A `-yes` never clears
  `plasti_off`.

## Verification (Professional binary 02-08-2026 as oracle)

quad4 tension_direct: `group_plasti_apply 0 -no` sigyy 21.978022 / eppyy 0;
`group_plasti_apply 0 -yes` + `plasti_apply -no` 21.978022 / 0. 1D bar
tension: 3 / 0; with visco_power: 3 / 0. Two bars, group 1 -no: 1 and 3.
GNU equal in every case. Tests: `group_plasti_apply1`, `group_plasti_apply2`.

## Pending

- `slide.cc` slide springs (gated by `plasti_apply`) and the
  truss/spring/interface plasticity are not gated per group (not measured).
