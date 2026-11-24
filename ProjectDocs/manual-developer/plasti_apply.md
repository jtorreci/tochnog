# plasti_apply / control_plasti_apply — developer

## Implementation

- **Enums**: `CONTROL_PLASTI_APPLY`, `PLASTI_APPLY` in `tochnog.h` /
  `tochnog-mod.h` (in sync), after `MATERI_FAILURE_APPLY` in the CONTROL
  block.
- **Records** (`database.cc`): `control_plasti_apply` INTEGER 1, class
  CONTROL (indexed); `plasti_apply` INTEGER 1, class CONTROL, `no_index 1`.
- **Gate**: `materi_apply_gate_off( PLASTI_APPLY, CONTROL_PLASTI_APPLY )`
  (`general.cc`, shared with `materi_damage_apply`): control record at the
  current `ICONTROL` overrules the global record; default `-yes`.
- **Consumers**:
  - `stress.cc`, `set_stress()`: `plasti_off` is evaluated once per call.
    1. incremental models (`plasti_incr`: tskh, aitskh,
       camclay_incremental): `plasti_incremental` forced false, so the
       classical loop runs with the elastic `C`;
    2. hypoplasticity: the `hypoplasticity()` call in the membrane loop is
       skipped (stress = elastic part only, normally none -> the rotated old
       stress; history untouched);
    3. direct cut-offs: `gate_off` of the pressure/coord-limit block is
       initialised with `plasti_off`, so `materi_direct_full_mc`,
       `materi_direct_cutoff` and `materi_compression_cutoff` are skipped;
    4. yield test: before `plasti_rule( ..., GET_YIELD_RULE, ... )` the loop
       sets `plasti_found = 1; continue;` (elastic trial state accepted,
       `inc_epp` = 0). All laws in `plasti_rule` (plasti.cc: mohr_coul,
       vonmises, druck_prag, tension, compression, camclay, cap1, cap2,
       hardsoil, diprisco, gurson, ...) are therefore off, and history
       variables written by `plasti_rule` keep their current values.
  - `slide.cc`, `slide_spring()`: without `control_slide_plasti_apply` the
    plastic cap follows `plasti_apply` / `control_plasti_apply` (manual
    6.371 "Default ... plasti_apply or control_plasti_apply will be used").
- Old plastic strains (`old_epp`) are untouched: a later `-no` timestep
  continues elastically from them (matches the Professional, see below).

## Verification (Professional binary 02-08-2026 as oracle)

No runnable Professional test covers the switches (`test/other/safety1/2.dat`
use `control_plasti_apply` inside the unimplemented safety_slip feature,
`validation_14_mesh.dat` is an include fragment). Own model: one plane-strain
quad4, E 1000, nu 0.3, `-total_linear`, top edge compressed with vely -0.02
(or +0.02 for tension laws) for 1 s, 10 steps; values at node 4.

| Law | Pro (a) data | Pro (b) `-no` | Pro (c) removed | GNU (a) | GNU (b) | GNU (c) |
|-----|--------------|---------------|-----------------|---------|---------|---------|
| mohr_coul 0.5 1 0 | sigyy -3.39353 | -21.978 | -21.978 | -3.48738 | -21.978 | -21.978 |
| vonmises 5 | -5.77357 | -21.978 | -21.978 | -5.77447 | -21.978 | -21.978 |
| druck_prag 0.5 1 0 | -6.87463 | -21.978 | -21.978 | -6.8168 | -21.978 | -21.978 |
| mohr_coul_direct 0.5 1 0 | -2.84076 | -21.978 | -21.978 | -2.43657 | -21.978 | -21.978 |
| tension_direct 1 (tension) | 1 | 21.978 | 21.978 | 1 | 21.978 | 21.978 |
| tension 1 (tension) | 1 | 21.978 | 21.978 | 1 | 21.978 | 21.978 |

Professional distribution tests with `plasti_apply -no` appended:

| Test | Quantity | Pro active | Pro `-no` | GNU active | GNU `-no` |
|------|----------|-----------|-----------|------------|-----------|
| other/hypo1.dat (wolffersdorff) | sigyy / hyhis0 | -862.929 / 0.598063 | -200 / 0.6 | -862.927 / 0.598045 | -200 / 0.6 |
| same, node 4 | disx | 0.0387914 | 2.3e-11 | — | 0 |
| validation/validation_13.dat (camclay) | cchis1 | 262.32 | 207 | 263.847 | 207 |
| other/hardsoil1.dat | sigyy node 1 | -369.013 | -7252.56 | GNU needs materi_plasti_kappa (not compared) | — |
| other/slide1.dat | velx node 10 | -0.00577061 | -0.0449775 (= control_slide_plasti_apply -no) | -0.00576771 | -0.0361781 (= control_slide_plasti_apply -no) |

Per-timestep (mohr_coul, two timesteps 0..0.5 / 0.5..1, sigyy):
control 1 `-no` Pro -14.3937 / GNU -14.4672; control 0 `-no` Pro -3.01552 /
GNU -2.46763; both `-no` -21.978 in both; `plasti_apply -no` +
`control_plasti_apply 1 -yes` = control 0 `-no`, `plasti_apply -yes` +
control 1 `-no` = control 1 `-no` in both binaries (precedence: control
overrules global). The residual Mohr-Coulomb differences are the
pre-existing return-mapping differences of case (a), not the gate. With
the exact tension_direct cut-off both binaries agree to all digits
(control 1 `-no`: 11.989; control 0 `-no`: 1; global `-no` + control 0
`-yes`: 11.989; global `-no` + control 1 `-yes`: 1).

Self-checking inputs (`validation-suite/test-2026/`): `plasti_apply1.dat`,
`plasti_apply2.dat`, `plasti_apply3.dat`, `control_plasti_apply1.dat`,
`control_plasti_apply2.dat`; rc=0 on both binaries, a changed target makes
the GNU return rc=1.

## Left out / pending

- `group_truss_plasti` (truss.cc), `group_spring_plasti` (spring.cc),
  `group_interface_materi_plasti_*` (interface.cc) and beam plasticity are
  not gated: their behaviour under `plasti_apply -no` was not measured on
  the Professional.
- `user_sigma` / `stress_umat` (user laws) are not gated.
- `group_plasti_apply` (6.762, per element group) and
  `linear_calculation_apply` (6.797, which sets `plasti_apply -no` among
  other things) are not implemented.
- `group_materi_plasti_hardsoil` c*cot(phi) is still used by the hardsoil
  elastic stiffness with `-no` (cohesionless in the only Pro example, so
  not discriminated).

## Differences with the Professional version

- **plasticity return mapping (MC, DP + refine, von Mises + refine)** (known
  difference, not decided; AUDIT K15). GNU: `validation_14` velx 48 % off;
  `validation_3` kap 11 % off; `mohrcou2` -6.932; MC `plasti_apply` model
  sigyy -3.487. Professional 02-08-2026: `mohrcou2` 0; MC sigyy -3.394.
  Status: UNKNOWN root cause; investigate before touching code. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`linear_calculation_apply`; truss/spring/interface/beam plasticity
  gating** (known difference, not decided; manual pages only). GNU:
  `linear_calculation_apply` not implemented; `plasti_apply` does not gate
  truss/spring/interface/beam plasticity. Professional 02-08-2026:
  `linear_calculation_apply` sets `plasti_apply -no`; gating of the other
  elements not measured. Status: Not implemented / not measured. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
