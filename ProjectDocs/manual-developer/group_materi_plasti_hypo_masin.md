# group_materi_plasti_hypo_masin

## Files and functions

- `masin.c` — **new pure-C port** of the reference Fortran `umat_hcea.for`
  (Tamagnini, Sellari, Masin, von Wolffersdorff; GPL). ~1060 lines.
  - `masin_umat()` — single integration-point update (Abaqus UMAT interface):
    `stress[6]`, `statev[16]`, `ddsdde[36]`, `dstran[6]`, `dtime`, `props[29]`,
    `testing`, `error`.
  - Internal: `get_tan()` (constitutive tangent L/H/M/N), `get_F_sig_q()`,
    `rhs()`, `rkf23_update()` (adaptive Runge-Kutta-Fehlberg 2/3 with
    substepping), `check_RKF()` (admissibility), `perturbate()` (tangent),
    `norm_res()`, `inv_sig()`, `inv_eps()`, `solout()`, `calc_statev()`.
  - Validated against the Fortran reference to `5e-7` on three strain paths
    (isotropic, triaxial, anisotropic; see `validation-suite/reference-masin/`).
- `hypoplas.cc` — dispatch block for `GROUP_MATERI_PLASTI_HYPO_MASIN`:
  reads parameters, converts tensors 3x3(row-major) to/from Voigt6, maps
  `hisv` to/from `statev`, calls `masin_umat()`, writes back
  stress/tangent/history.
- `tochnog.h` / `tochnog-mod.h` — enum entries:
  `GROUP_MATERI_PLASTI_HYPO_MASIN`, `_STRUCTURE`, `_OCR`,
  `CONTROL_MATERI_PLASTI_HYPO_MASIN_OCR_APPLY` (same order, must stay in sync).
- `database.cc` — keyword registrations:
  - `group_materi_plasti_hypo_masin`: DOUBLE, length 5.
  - `group_materi_plasti_hypo_masin_structure`: DOUBLE, length 3.
  - `group_materi_plasti_hypo_masin_ocr`: DOUBLE, length 1.
  - `control_materi_plasti_hypo_masin_ocr_apply`: INTEGER, CONTROL class.
- `check.cc` — requires `materi_stress` and `materi_history_variables`.
- `Makefile` — `MASIN_SRC=masin.c`, `MASIN_OBJ=masin.o`, compiled like
  `hypo.c` (pure C, no f2c); `-lm` already in the link line.

## Implementation details

- **Parameter mapping** (the 5-record values do NOT map 1:1 to `props`):
  `group_materi_plasti_hypo_masin` = `phi_c lambda* kappa* N r` is read into a
  temp array, then `props[0]=phi_c, props[2]=lambda*, props[3]=kappa*,
  props[4]=N, props[5]=r`. `props[1]=p_t` is reserved for the cohesion shift
  (kept 0). `GET_AND_CHECK` needs the explicit length (5).
- **Defaults** applied: `props[6]=1` (alpha_G isotropic), `props[9]=1`
  (s_f), `props[13]=0` (A_g: intergranular strain off), `props[17]=3`
  (vertical direction z). Structure record fills `props[7,8,9]=k,A,s_f`.
- **History layout** (`materi_history_variables >= 8`):
  `hisv[0..5]`=intergranular strain, `hisv[6]`=e, `hisv[7]`=sensitivity.
  The Fortran reference defines `move_asv_hcea` (which negates the
  intergranular strain) but NEVER calls it in the integration path
  (`iniy_hcea` copies `asv` directly); therefore no sign flip is applied.
- **Tensor conversion** (3x3 row-major -> Voigt6):
  `mstress[0]=sig[0], [1]=sig[4], [2]=sig[8], [3]=sig[1], [4]=sig[2],
  [5]=sig[5]`; symmetric on output. Same for strain.
- **Tangent conversion** (Voigt6 ddsdde -> 3x3x3x3 Chypo): uses the Voigt
  index map `vi[0][0]=0, vi[1][1]=1, vi[2][2]=2, vi[0][1]=3, vi[0][2]=4,
  vi[1][2]=5` (symmetric pairs).
- **Stress accumulation** follows the `hypo_` pattern:
  `new_sig += stress_total - rotated_old_sig` (incremental formulation;
  `formulation==TOTAL` is rejected).
- The kernel is called with `testing=0` (strict tolerance `tolintT=1e-3`,
  `maxnint=10000`). The internal adaptive substepping resolves the strain
  increment regardless of the global `control_timestep`.

## Validation

- Kernel: `validation-suite/reference-masin/` holds the Fortran reference
  driver (`masin_ref_driver.f90`), the UMAT source (`umat_hcea.for`) and the
  generated reference strain paths (`iso_path1.txt`, `iso_path2.txt`,
  `aniso_path2.txt`). The C port reproduces them to `5e-7` (the residual is
  only the Fortran print precision `E20.10`).
- End-to-end: `validation-suite/test-2014/hypomasin1.dat` (single quad4,
  laterally confined biaxial path). Targets: `sigxx=-340.7` (reference
  -334.8, +1.8%), `hisv6=0.6333` (reference 0.6663, -5%). The difference
  comes from the strain path of the single-element FE setup vs the ideal
  driver path, NOT from the kernel.
- Regression: hypo1-4 still pass (wolfersdorff unaffected).

## External dependencies

- `masin.c` is self-contained (only `<math.h>`, `<stdio.h>`, `<stdlib.h>`);
  no f2c, no tochnog headers. The `masin_umat` prototype is declared
  `extern "C"` in `hypoplas.cc`.
- Core `db()`, `db_active_index()`, `pri()`, `array_add`/`array_subtract`
  from tochnog.

## Hardcoded parameters / pending refactorings

- `mtesting=0` is hardcoded; the PLAXIS first-call loose-tolerance mode
  (`testing=1`) and the stiffness-only mode (`testing=2`) are not exposed.
- The intergranular strain slots (`hisv[0..5]`) are reserved but inactive
  (`A_g=0` forces the basic law); activating them requires `props[13]>0`
  plus reading `group_materi_plasti_hypo_strain_intergranular_masin_clay`.
- The clay anisotropic variant (alpha_G, alpha_E, alpha_nu, direction) is
  hardcoded to isotropic defaults; the advanced/clay keywords are not
  registered yet (P4-B2).
- `materi_history_variables >= 8` is enforced with an explicit error; the
  check could be lifted to a softer warning for backwards compatibility.
- The `OCR` initial-void-ratio formula duplicates the Fortran; verify against
  the Fortran when the OCR path is exercised.
- `masin.c` uses fixed `props[29]` and `statev[16]` sizes; refactor to
  structs if more variants (strength reduction, visco) are ported.
