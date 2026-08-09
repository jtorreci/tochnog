# group_materi_hyper_besseling (and the hyperelastic family)

## Files and functions

- `hyperela.cc` — the whole hyperelastic family lives here:
  - `hyperelasticity()` (line 26) — entry point: computes `sig` from `epe` and,
    when `group_materi_hyper_stiffness -yes` (default), also fills `Chyper`.
  - `hyper_Cmat()` (line 46) — consistent tangent by central differences of the
    stress with respect to the strain (`DELTA=1.e-2`, `EPS_VARIATION=1.e-3`).
    Diagonal entries `tmp` are clamped to `>=0`.
  - `hyper_stress()` (line 82) — builds `C = F^T.F` from `epe` depending on
    `memory`:
    - `-TOTAL_PIOLA`: `C = 2*epe + I`
    - `-TOTAL`, `-TOTAL_LINEAR`, or when `materi_strain_elasti` is active:
      `C = (epe + I)^T . (epe + I)`
    - otherwise returns 0 (no stress).
    It then computes the stress by central differences of the strain energy `W`.
  - `hyper_law()` (line 137) — evaluates `W` for every active
    `GROUP_MATERI_HYPER_*` record and sums the contributions. `W<0` is clamped
    to 0 (line 225).
- Entry point called from `stress.cc:664-667` (both old and new configuration)
  inside the `materi_stress` block, under the `hyperelasticity` branch.
- `database.cc:2388-2443` — keyword registrations for all 10 hyper records plus
  `GROUP_MATERI_HYPER_STIFFNESS` (line 2419, INTEGER, default `-YES` in
  `hyperela.cc:31`).
- Enums mirrored in `tochnog.h` (lines ~539-548) and `tochnog-mod.h`.

## Implementation details

- All models are evaluated through ONE strain-energy function; the active
  models are detected with `get_group_data(..., GET_IF_EXISTS)` inside
  `hyper_law()`, so several hyper records can be combined (energies add).
- `db_error()` is raised for invalid parameter combinations: `beta=0`/`G=0`
  (Blatz-Ko), `beta==1` (Murnaghan), `beta==0` (Ogden). A pathological `J==0`
  aborts with "Probably too large element distortion".
- The stress is NOT an explicit derivative of `W`; it is always computed by
  central differences (`2*(W_right-W_left)/(2*variation)`), which is why the
  consistent tangent `Chyper` is also obtained by differences of `stress` — the
  whole family is purely numerical at the constitutive level.
- `group_materi_hyper_stiffness -no` disables the `Chyper` computation but
  keeps the stress; used to save CPU time when the elastic tangent is an
  acceptable approximation.
- Regression coverage: `validation-suite/test-2014/blatz1.dat` (Blatz-Ko,
  target `sigyy=3.18021±0.01` at stretch 2.0), `ho_mech1.dat`
  (Mooney-Rivlin, combined with a plasticity model).

## External dependencies

- Core `db()`/`get_group_data()`/`db_error()` accessors only.
- `matrix_invariants()` (miscel.cc) — I1, I2, I3 of `C`.
- Globals `materi_strain_elasti`, `MDIM`, `TOTAL_PIOLA`, `TOTAL`,
  `TOTAL_LINEAR` from `tochnog.h`.

## Hardcoded parameters / pending refactorings

- `DELTA=1.e-2` and `EPS_VARIATION=1.e-3` are hardcoded (hyperela.cc:23-24):
  they control the finite-difference steps for both stress and tangent. They
  are a compromise between numerical noise and truncation error; there is no
  user-facing record to tune them.
- The tangent diagonal clamp (`if tmp<0 tmp=0`) is a crude
  positive-definiteness fix and may mask real softening in the `kdim==ldim`
  entries.
- `W` is clamped to `>=0`; this makes the tangent zero under severe
  compression where the exact energy would turn negative.
- The enum + `strcpy` registration is duplicated between `tochnog.h` and
  `tochnog-mod.h` (project-wide convention; a single shared enum would remove
  the drift risk).
