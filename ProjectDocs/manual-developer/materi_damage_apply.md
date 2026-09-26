# materi_damage_apply — developer

## Implementation

- **Enum**: `MATERI_DAMAGE_APPLY` in `tochnog.h` / `tochnog-mod.h` (in sync),
  right after `CONTROL_MATERI_DAMAGE_APPLY` (same placement pattern as
  `MATERI_DYNAMIC` after `CONTROL_MATERI_DYNAMIC`).
- **Record** (`database.cc`, `db_initialize`): INTEGER, length 1,
  `no_index 1`, data_class CONTROL.
- **Helper** (`general.cc`): `materi_apply_gate_off( global_item,
  control_item )` returns 1 when the feature must be ignored for the current
  timestep. It reads the global record (index 0, `GET_IF_EXISTS`), then the
  control record at the current `ICONTROL` (`GET_IF_EXISTS`), which
  overrules it; default `-yes`. Declared in both headers next to
  `control_materi_gate_off` (the older control-only helper, kept for the
  other `control_materi_*_apply` gates).
- **Consumer** (`damage.cc`, `damage()`): `materi_apply_gate_off(
  MATERI_DAMAGE_APPLY, CONTROL_MATERI_DAMAGE_APPLY )` -> `new_damage =
  old_damage` and return (damage frozen, not reset). The caller in
  `set_stress` (stress.cc) still scales the stress by `1-new_damage`, so a
  damage reached in an earlier timestep keeps acting.
- **Behaviour change**: the previous gate also tested
  `CONTROL_MATERI_FAILURE_APPLY`; it was removed because the Professional
  02-08-2026 does not stop damage with the failure switch (see
  `materi_failure_apply.md`).

## Verification (Professional binary 02-08-2026 as oracle)

Model: the Professional distribution test `test/other/damage1.dat` (1D bar,
Mazars), targets `dam` and `sigxx` at node 2.

| Case | Pro dam | Pro sigxx | GNU dam | GNU sigxx |
|------|---------|-----------|---------|-----------|
| (a) damage active | 0.993903 | 85598.3 | 0.993903 | 1436.83 (pre-existing difference) |
| (b) `materi_damage_apply -no` | 0 | 1.404e7 | 0 | 1.404e7 |
| (c) damage record removed | 0 | 1.404e7 | 0 | 1.404e7 |
| `control_materi_damage_apply 0 -no` | 0 | 1.404e7 | 0 | 1.404e7 |
| global `-no` + control 0 `-yes` | 0.993903 | 85598.3 | 0.993903 | 1436.83 |
| global `-yes` + control 0 `-no` | 0 | 1.404e7 | 0 | 1.404e7 |
| `materi_failure_apply -no` | 0.993903 | 85598.3 | 0.993903 | 1436.83 |

Two timesteps (0..3e-4, 3e-4..6e-4): control 1 `-no` -> dam 0.451188 (frozen)
in both; control 0 `-no` -> 0.993903 in both; global `-no` + control 1 `-yes`
-> 0.993903; global `-no` + control 0 `-yes` -> 0.451188 in both.

Precedence: control record of the current index overrules the global one.

Self-checking inputs `validation-suite/test-2026/materi_damage_apply1.dat`
and `materi_damage_apply2.dat` pass on both binaries; changing a target
value makes the GNU return rc=1.

## Pending

- The damaged stress (`(1-d)` scaling path) differs from the Professional
  independently of this switch (e.g. at 3e-4: GNU sigxx 1.68071e6, Pro
  3.85266e6 with the same dam 0.451188); not addressed here.
