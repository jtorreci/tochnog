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
- **materi.cc, `materi()`** — rotation to the new configuration: included
  in the `inc_rot` branch of `-UPDATED`/`-UPDATED_LINEAR` (rot is the
  identity because of the linear branch below, so no rotation is applied).
- **materi.cc, `set_deften_etc()`** — identity rotation matrices and
  LINEAR engineering strains: `-UPDATED_AREA` added to both conditions.
- **stress.cc, `set_stress()`** — compressibility block: `-UPDATED_AREA`
  added to the linear-strain formula (`tmp = tr(epe)/compressibility`),
  same as `-TOTAL_LINEAR`/`-UPDATED_LINEAR`.

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
Until it exists, `-updated_area` and `-updated_linear` are behaviorally
identical in the GNU. The four corpus tests using `-updated_area` are
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
