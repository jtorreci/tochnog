# group_materi_elasti_k0 (+ control_materi_elasti_k0 hook)

## Implementation

- **Enum**: `GROUP_MATERI_ELASTI_K0` in tochnog.h and tochnog-mod.h (in
  sync, alphabetical between `ELASTI_COMPRESSIBILITY` and `ELASTI_LADE`).
- **Record** (database.cc): DOUBLE_PRECISION, length 1, data_class MATERI,
  data_required GROUP_TYPE: `K0`. `CONTROL_MATERI_ELASTI_K0` (INTEGER,
  length 1, CONTROL) is the per-timestep-index switch.
- **Hook in `set_stress()`** (stress.cc), before the young block:
  ```c
  k0_active = 0;
  k0_control_swit = 0;
  idum[0] = 0;
  db( ICONTROL, 0, idum, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_MATERI_ELASTI_K0, idum[0], &k0_control_swit, ddum, ldum,
    VERSION_NORMAL, GET_IF_EXISTS );
  if ( k0_control_swit==-YES &&
       get_group_data( GROUP_MATERI_ELASTI_K0, gr, element, new_unknowns,
         &k0_elasti, ldum, GET_IF_EXISTS ) ) {
    if ( k0_elasti>0.95 ) k0_elasti = 0.95;
    k0_active = 1;
  }
  ```
- **young / young_power blocks**: the poisson record is read (if present)
  and then `if ( k0_active ) poisson = k0_elasti/(1.+k0_elasti);` — the
  override no longer depends on the poisson record being present.
- **hardsoil block**: right after reading the 7 hardsoil parameters,
  `if ( k0_active )` sets `nu50 = nuur = K0/(1+K0)`, so the trial C of the
  loading/unloading decision and the final C both use the K0 Poisson.
- **check.cc**: `GROUP_MATERI_ELASTI_K0` requires `materi_stress` +
  `materi_velocity`.

## History of fixes (2026-09-29)

1. **Index guard**: the hook was wrapped in
   `if ( db_active_index( CONTROL_MATERI_ELASTI_K0, 0, VERSION_NORMAL ) )`.
   `db_active_index` tests that EXACT index, so a switch given for any
   control index other than 0 (Professional test `elasti7`: index 20) was
   silently ignored (sigxx 0 instead of -0.05). The record is now read
   directly at the ICONTROL index with GET_IF_EXISTS (`k0_control_swit`
   reset to 0 because GET_IF_EXISTS does not write when absent). Inputs
   with index 0 are byte-identical.
2. **Poisson record not required**: measured on Professional 02-08-2026:
   young + K0 without `group_materi_elasti_poisson`, and young_power + K0
   with and without poisson, all give sigxx = K0*sigyy = -0.05.
3. **hardsoil**: the old comment claimed hardsoil was NOT affected; the
   manual (6.650) lists it and the Professional confirms it: hardsoil
   (`nu50 = nuur = 0` or `0.2`) + K0 0.1 gives sigxx -0.05, and without the
   control switch 0 / -0.125.
4. K0 > 0.95 clamp verified: K0 = 2 gives sigxx/sigyy = 0.95 on both
   binaries.

## Validation

| Test | Before | After | Professional |
|------|--------|-------|--------------|
| `elasti7` (index 20) sigxx / sigyy / top disy | 0 / -0.5 / -5e-7 | -0.05 / -0.5 / -4.9090909e-7 | -0.04999999884 / -0.5 / -4.909090913e-7 |
| `elasti_k0_index1` (test-2026) | rc 1 (sigxx 0) | PASS | -0.05 |
| `elasti_k0_nopoisson1` (young_power, no poisson) | rc 1 (sigxx 0) | PASS | -0.05 |
| `elasti_k0_hardsoil1` (hardsoil nu 0.2) | rc 1 (sigxx -0.125) | PASS | -0.05 |
| `mk0` / `mk0_off` (index 0) | -0.75 / -0.5769 | unchanged | — |

All elasti7 variables (node_dof, post_point_dof, node_rhside) agree with
the Professional within 2.4e-8 relative after the fix.

## Pending / notes

- young_polynomial does not receive the override (the manual only lists
  young/young_power/hardsoil) — not measured on the Professional.
- The GNU input check requires `materi_strain_total` for
  `group_materi_elasti_young_power`; the Professional runs without it
  (separate legacy check, not changed here).
