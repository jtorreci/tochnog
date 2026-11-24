# group_materi_memory -updated_area

## Semantics (Professional)

`-updated_area` (manual Professional 6.784, incremental_driver section):
small deformation theory (same kinematics as `-updated_linear`) with the
area change of the loaded surface taken into account. It is a driver
option: it exists to make force-controlled laboratory experiments
(oedometer / triaxial / direct shear) apply the prescribed stress on the
CURRENT cross-section of the deforming specimen.

## Implementation

The option is a new value of the material memory enum:

1. `UPDATED_AREA` added to the enum of `tochnog.h` and `tochnog-mod.h`
   (same position in both: between `UPDATED` and `UPDATED_LINEAR`).
2. `name[UPDATED_AREA] = "updated_area"` registered in the name table of
   `database.cc` (next to `updated`/`updated_linear`), so the value
   parses in any memory record (`group_materi_memory` and the structural
   `group_*_memory` records).

The kinematic treatment follows `-updated_linear` branch by branch:

- **materi.cc, `materi()`** — back-rotation of the old stress: included in
  the direct `array_move` branch (memory==-UPDATED ... || -UPDATED_AREA).
- **materi.cc, `materi()`** — rotation to the new configuration:
  `-UPDATED_AREA` is NOT in the rotation list (only `-UPDATED`, `-TOTAL`,
  `-TOTAL_PIOLA` rotate). CORRECTION 2026-09-29: this page used to say that
  `-UPDATED_AREA` sat in the `inc_rot` branch and that `inc_rot` was the
  identity. That was wrong: `set_deften_etc()` sets `inc_rot` to the
  identity for the linear memories but then calls
  `set_deften_u_rot(inc_deften, inc_u, inc_rot)` unconditionally, which
  overwrites it with the real incremental rotation of the polar
  decomposition. The stress was therefore rotated by the incremental
  rotation (simple shear, gamma = 1: rotated by atan(gamma/2) = 26.57 deg).
  Removed from the list together with `-UPDATED_LINEAR` (audit cluster K01,
  see `group_materi_memory_updated_linear.md`).
- **materi.cc, `set_deften_etc()`** — identity `old_rot`/`new_rot` and
  LINEAR engineering strains: `-UPDATED_AREA` added to both conditions.
  The velocity gradient is the CURRENT-geometry one (the reference-geometry
  gradient of commit `c945833` is gated on `-UPDATED_LINEAR` only), so on a
  `-follow_material` mesh the accumulated strain is the logarithmic-like
  measure, not the engineering one (1D bar compressed 0.5%: `-updated_area`
  -0.0050117, `-updated_linear` -0.005).
- **stress.cc, `set_stress()`** — compressibility block: `-UPDATED_AREA`
  takes the linear-strain branch, same as `-TOTAL_LINEAR`/`-UPDATED_LINEAR`.
  Since 2026-09-29 the incremental volumetric stress is
  `tr(inc_epe)/compressibility` (it was `tr(new_epe - old_epe)`, which adds
  the total strain every step when `materi_strain_elasti` is not
  initialized).

NOT touched on purpose:

- **stress.cc `control_materi_updated_apply` -no conversion**: converts
  only the large-deformation memories (`-UPDATED`, `-UPDATED_WITHOUT_ROTATION`)
  to `-updated_linear`; `-updated_area` is already small deformation, so no
  conversion applies (no corpus test combines both records; revise when the
  incremental_driver lands).
- **stress.cc displacement check** (memory==-UPDATED ||
  -UPDATED_WITHOUT_ROTATION rejects `materi_displacement`): `-updated_area`
  is not included, mirroring `-updated_linear`.
- **elem.cc new_coord update**: coordinates are updated for every memory
  except `-TOTAL_LINEAR`, so `-updated_area` integrates on the current
  geometry (this is the side of the semantics the GNU already covers).
- **Structural records** (`group_beam_memory`, `group_contactspring_memory`,
  `group_truss_memory`, `group_interface_materi_memory`): their value
  lists are unchanged; `-updated_area` on those records falls into their
  existing default/`db_error` handling (not used by the corpus).

## Pending

The load-area part of the semantics (stress x CURRENT area in
force-controlled experiments) has NO carrier in the GNU: the
`incremental_driver` machinery (section 6.784: experiments, equilibrium
loops, force/displacement control, `incremental_driver_result.txt`
output, groundflow coupling for undrained tests) is not implemented.
Until it exists, `-updated_area` and `-updated_linear` differ in the GNU
only by the velocity gradient (current geometry for `-updated_area`,
reference geometry for `-updated_linear`, see above); rotation handling,
strain measure and compressibility are the same.

The Professional 02-08-2026 REJECTS `-updated_area` in an ordinary
`group_materi_memory` record ("Error detected for data item :
group_materi_memory"), measured 2026-09-29 on the sfnet `druckpr1` and Pro
`interac1` models with the memory switched to `-updated_area`. Its manual
(6.685) lists only `-updated`, `-updated_jaumann`, `-updated_linear`,
`-total`, `-total_linear`; `-updated_area` appears only in the
`incremental_driver` section. So the GNU-only use of `-updated_area` outside
the incremental driver has no Professional oracle. The four corpus tests using `-updated_area` are
blocked by the `incremental_driver` family, NOT by this option (verified:
with the value accepted, the parse of the test files proceeds until the
first unknown record `incremental_driver`).

## Verification

- Build clean + own suite 16/16.
- A/B on the GNU validation model hypo2 (hypoplasticity + intergranular
  strain, oedometric compression): `-updated_area` vs `-updated_linear`
  both rc=0, identical output and identical check targets (sigxx
  -0.17361, hisv0 0.63843).
- Corpus parse progression: the four `incremental_driver/` tests now fail
  at the first unknown record `incremental_driver` instead of the value
  `-updated_area` (oedometric/triaxial drained/undrained: RUNFAIL;
  incremental_driver_syntax: PARSE, blocked earlier in initia by
  `materi_strain_isa_c`).

### 2026-09-29 (audit cluster K01 fix)

GNU after the fix, `-updated_area` probes (no Professional oracle, see
above): sfnet `druckpr1` (simple shear, DP c = 0.1) sigxx 0, sigxy 0.11547,
sigyy 0 (pre-fix 0.0924 / 0.0693 / -0.0924); Pro `interac1` (2-bar,
compressibility) sigxx -0.0050117 at both ends (pre-fix -0.010047, element 2
strain ~0). The residual 0.23% against the analytic -0.005 is the
current-geometry gradient on the GNU's default `-follow_material` mesh.
