# control_materi_updated_apply

## Implementation

- Record: `database.cc` (`CONTROL_MATERI_UPDATED_APPLY`, INTEGER, class
  CONTROL, length 1).
- `-yes` direction: `general.cc` `materi_memory_default()` returns
  `-UPDATED` instead of `-UPDATED_LINEAR` when the record at the current
  `ICONTROL` is `-YES`. Only an unspecified memory is affected, because
  every reader assigns the default before reading the explicit
  `GROUP_MATERI_MEMORY` record (materi.cc, stress.cc, elem.cc, tendon.cc).
- `-no` direction: `stress.cc` `set_stress()`, after the memory read:
  `control_materi_gate_off(CONTROL_MATERI_UPDATED_APPLY)` turns
  `-UPDATED` / `-UPDATED_WITHOUT_ROTATION` into `-UPDATED_LINEAR`
  (Sprint 9; only in set_stress, not in the kinematics of elem/materi).
- Index semantics: the record with the index of the current timestep
  (`ICONTROL`), as the other control_materi_* gates.

## Verification (2026-10-04)

Single quad4 simple shear gamma 1, velocity kinematics, no memory record:

| switch | GNU | Professional 02-08-2026 |
|---|---|---|
| none | 0 / 0.5 / 0 | 0 / 0.5 / 0 |
| `0 -yes` | 0.231098 / 0.420739 / -0.228598 | 0 / 0.5 / 0 |
| explicit `-updated` + `0 -yes` | 0.231098 / 0.420739 / -0.228598 | same |
| explicit `-updated` + `0 -no` | 0 / 0.5 / 0 | 0 / 0.5 / 0 |

The Professional ignores `-yes` also with `materi_displacement`, an
explicit `mesh -follow_material` and `mesh -fixed_in_space`, and for
index 1. Decision 2026-10-04: follow manual 6.152 (deliberate
divergence). Test: `validation-suite/test-2026/materi_memory_default_updated_apply1.dat`
(the Professional fails its targets). No corpus input uses `-yes`.

## Pending

- The `-no` conversion is applied only in `set_stress()`; the element
  kinematics (elem.cc / set_deften_etc) still see the original memory.
