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
- `masin_visco.c` — **new pure-C port** of the visco extension
  `umat_visco.f` (Jerman & Masin 2020, GPL). Same skeleton as `masin.c` but:
  - `NASV=10` (2 extra state slots; `materi_history_variables >= 10`).
  - `inv_sig()` computes rotated invariants (`sig_rot` from the shear-band
    angle `beta`, `cos3t_rot`, `I1rot/I2rot/I3rot`).
  - `get_tan()` uses `Fmfactor_rot` (rotated Matsuoka-Nakai), `ocparam`
    (parms[21]) instead of ocrcs=2, and the LD approach: `gama`
    (flow-direction angle, auto-initialized from `tangama` when < -pi/2),
    `hypo_Dsom_ld`, `wy`/`acorrwy`, and the `LL_unl`/`LL` split.
  - `get_F_sig_q()` divides `deps` by `dtime` (rate formulation), scales by
    `Dref`, and re-multiplies — the rate-dependent visco response.
  - `masin_visco_umat()` reduces `phi_c` by `beta` (shear-band softening)
    and reads e0/OCR from `props[27]`.
  - Validated against `umat_visco.f`: identical `sig11=-82.5684` at step 20.
  - `masin_niemunis_visco_umat()` — the **Niemunis visco law** (`Dr Iv`)
    implemented from the professional manual theory (no reference UMAT
    exists; the soilmodels "Niemunis" page is the unrelated High-Cycle
    Accumulation model). Computes `L = fb*Lhat`, the flow rule `m`, the
    Niemunis `OCR = pe/pe+` and the creep rate `Dr*(1/OCR)^(1/Iv)`. The
    internal parameters default to clay values: `lambda=lambda*`,
    `ee0=e_initial`, `pe0=p_initial`, `betaR=1`. The creep exponent is
    clamped (`creep_rate <= 1e3*Dr`) to keep the simplified defaults stable.
- `masin_basic.c` — **new pure-C kernel** for the five-parameter Masin
  (2005) basic route, from Appendix C, Eqs. (73)-(85). Deliberately separate
  from `masin.c` (HCEA) and `masin_visco.c`.
  - `masin_basic_umat()` — Abaqus-style update: `stress[6]`,
    `statev[6]` = void ratio `e`, `statev[13]` = sensitivity `s`, `ddsdde[36]`,
    `dstran[6]`, `props[29]`. The kernel reads `props[0]=phi_c, [2]=lambda*,
    [3]=kappa*, [4]=N, [5]=r` and the structure parameters
    `props[7]=k, [8]=A, [9]=s_f`.
  - Internal: `mb_dot()`/`mb_det()` (Voigt6 metric with the factor 2 on the
    shear slots), `mb_rate()` (the Appendix C rate law plus the structure
    terms), `mb_integrate()` (fixed 32-substep third-order Runge-Kutta, now
    integrating `e` and `s`), and a numerical algorithmic tangent (central
    difference, `eps=1e-7`).
  - Source-grounded fixes: (1) the m-tensor Eq. (15) uses `T_hat:T_hat` in
    the numerator and denominator (`that2 = mb_dot(that, that)`), not the
    `T_hat*` variants; (2) the `cos3theta` invariant `tr(T_hat^3)`
    (Eq. 82 / manual Eq. 17) expands its cross terms with coefficient 3
    (`3*a12^2*(a11+a22) + 3*a13^2*(a11+a33) + 3*a23^2*(a22+a33) +
    6*a12*a13*a23`), not 2.
  - **Structure (sensitivity) extension** (manual "Masin law" Eqs. 12-21 +
    Masin 2007 CGJ 44(5) Eq. 8): `structure_on = (s >= 1.0)`. When active,
    `mb_rate()` computes `Si = (s - k(s - s_f))/s` (manual Eq. 13), uses
    `kappa* Si` in `alpha` (Eq. 18), divides `c1` by `Si` (Eq. 19), multiplies
    `fs` by `Si` (Eq. 12) and divides the `fd` bracket by `s`
    (`fd = (2p/(s p*_e))^alpha`, 2007 Eq. 8; the manual prints the reference
    `s=1` form). `mb_integrate()` carries `s` through the same RK3 stages with
    `sdot = -(k/lambda*)(s - s_f) sqrt((trD)^2 + (A/(1-A))(2/3)||devD||^2)`
    (manual Eq. 21). When `structure_on` is false the code path is the exact
    original expression (`si = 1`, reference `alpha`/`c1`/`fs`/`fd`,
    `sdot = 0`), so the non-structure route is bit-identical.
  - `hypoplas.cc` dispatches the basic route
    (`!hypo_masin_clay && !hypo_masin_visco && !hypo_masin_visco_jm`) to
    `masin_basic_umat()`; `makefile` adds `MASIN_BASIC_SRC=masin_basic.c`
    and `masin_basic.o` to the link.
- `hypoplas.cc` — dispatch block for `GROUP_MATERI_PLASTI_HYPO_MASIN`:
  reads parameters, converts tensors 3x3(row-major) to/from Voigt6, maps
  `hisv` to/from `statev`, calls `masin_umat()`, writes back
  stress/tangent/history.
- `tochnog.h` / `tochnog-mod.h` — enum entries:
  `GROUP_MATERI_PLASTI_HYPO_MASIN`, `_STRUCTURE`, `_OCR`,
  `CONTROL_MATERI_PLASTI_HYPO_MASIN_OCR_APPLY`, `_MASIN_CLAY`,
  `_MASIN_CLAY_ADVANCED_PARAMETERS`, `_MASIN_CLAY_AVANCED_DIRECTION`,
  `_MASIN_CLAY_OCR`, `_MASIN_CLAY_STRUCTURE`,
  `CONTROL_MATERI_PLASTI_HYPO_MASIN_CLAY_OCR_APPLY`,
  `_HYPO_STRAIN_INTERGRANULAR_MASIN_CLAY` (same order, must stay in sync).
- `database.cc` — keyword registrations:
  - `group_materi_plasti_hypo_masin`: DOUBLE, length 5.
  - `group_materi_plasti_hypo_masin_structure`: DOUBLE, length 3.
  - `group_materi_plasti_hypo_masin_ocr`: DOUBLE, length 1.
  - `control_materi_plasti_hypo_masin_ocr_apply`: INTEGER, CONTROL class.
  - `group_materi_plasti_hypo_masin_clay`: DOUBLE, length 5.
  - `group_materi_plasti_hypo_masin_clay_advanced_parameters`: DOUBLE, length 4.
  - `group_materi_plasti_hypo_masin_clay_avanced_direction`: INTEGER, length 1.
  - `group_materi_plasti_hypo_masin_clay_ocr`: DOUBLE, length 1.
  - `group_materi_plasti_hypo_masin_clay_structure`: DOUBLE, length 3.
  - `control_materi_plasti_hypo_masin_clay_ocr_apply`: INTEGER, CONTROL class.
  - `group_materi_plasti_hypo_strain_intergranular_masin_clay`: DOUBLE, length 7.
- `check.cc` — requires `materi_stress` and `materi_history_variables`.
- `Makefile` — `MASIN_SRC=masin.c`, `MASIN_OBJ=masin.o`, compiled like
  `hypo.c` (pure C, no f2c); `-lm` already in the link line.

## Implementation details

- **Parameter mapping** (the 5-record values do NOT map 1:1 to `props`):
  `group_materi_plasti_hypo_masin` / `_clay` = `phi_c lambda* kappa* N r/nu_pp`
  is read into a temp array, then `props[0]=phi_c, props[2]=lambda*,
  props[3]=kappa*, props[4]=N, props[5]=r|nu_pp`. `props[1]=p_t` is reserved
  for the cohesion shift (kept 0). `GET_AND_CHECK` needs the explicit
  length (5).
- **Clay anisotropic mapping** (P4-B2a):
  - `_clay_advanced_parameters` = `alpha_G alpha_f ay oc` is read into
    `props[6]`, `props[20]`, `props[22]`, `props[23]` (in that order).
  - `_clay_avanced_direction` `diri` (0/1/2) maps to `props[17]=diri+1`.
  - `_clay_structure` fills `props[7,8,9]=k,A,s_f`.
  - `alpha_E`/`alpha_nu` stay 0 -> the kernel auto-derives them
    (`alpha_E=alpha_G^1.25`, `alpha_nu=alpha_G`).
  - In `masin.c`, `ay`/`oc` are read from `parms[22]`/`parms[23]` with
    defaults 0.30/2.0 when absent (both in `get_tan` and `check_RKF`).
- **Intergranular strain masin clay** (P4-B2b):
  `_strain_intergranular_masin_clay` = `R A_g n_g m_rat beta_r chi [theta]`
  maps to `props[10]=R, [13]=A_g, [14]=n_g, [15]=m_rat, [11]=beta_r,
  [12]=chi`. `A_g>0` activates the kernel `istrain=1` branch (small-strain
  stiffness). `theta` has no direct slot (the kernel interpolates with chi).
- **Defaults** applied: `props[6]=1` (alpha_G isotropic), `props[9]=1`
  (s_f), `props[13]=0` (A_g: intergranular strain off), `props[17]=3`
  (vertical direction z), `props[22]=props[23]=0` (ay/oc kernel defaults).
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
- End-to-end (single quad4, laterally confined biaxial path):
  - `hypomasin1.dat` (basic law): `sigxx=-340.7` (reference -334.8, +1.8%),
    `hisv6=0.6333` (reference 0.6663, -5%).
  - `hypomasin2.dat` (clay anisotropic, `alpha_G=2`): `sigxx=-429.6`
    (reference -418.6, +2.6%). The void ratio is identical to the isotropic
    case (anisotropy only affects the deviatoric response).
  - `hypomasin3.dat` (intergranular strain): `sigxx=-155.8`, `hisv6=0.6932`.
    NOTE: the kernel matches the Fortran reference (sig11=-195.4 at 20
    steps), but the FE end-to-end result differs (~20%) because the
    intergranular stiffness is sensitive to the equilibrium iterations
    inside each substep (each Newton iteration re-integrates with the delta
    of the previous iteration but the stress of the start of the step).
    The basic/anisotropic paths agree closely because they are not
    path-history-sensitive in the same way.
  - `hypomasin4.dat` (visco JM, `Dref=0.5`): `sigxx=-82.88` (Fortran
    reference -82.57, +0.4%), `hisv6=0.6932`. The JM kernel relaxes the
    stress under constant strain rate, matching the Fortran `umat_visco.f`
    (identical at the driver level: sig11=-82.5684).
  - `hypomasin5.dat` (Niemunis visco, `Dr=1e-6 Iv=0.1`): `sigxx=-104.1`.
    No numeric reference exists (theory-only implementation); the test pins
    the value and guards against regressions of the creep formulation.
- Regression: hypo1-4 still pass (wolfersdorff unaffected).
- **Basic route (Masin 2005), integrated with the R2 reference-geometry
  strain** (2026-09-13): end-to-end `hypo7`/`hypo9` (drained axisymmetric
  triaxial, `-updated_linear`) reproduce the Professional within the target
  tolerance. `sigyy` / `e`: hypo7 `-231.789` / `1.110832` vs Professional
  `-231.809` / `1.110746`; hypo9 `-231.789` / `1.110832` vs `-231.811` /
  `1.110740`; target `±0.1` → both PASS. The `ts3` cross-term fix is inert
  for the axisymmetric tests (zero shear).
- **Structure (sensitivity) evolution Eq. (21)** (`5a98958`, 2026-09-13):
  `hypo8` (Pisa clay, structure `k=0.4 A=0.1 s_f=1`, `e0=1.03`, `s0=2`)
  now passes: `rc=0`, `sigyy=-230.534`, `e=1.284535`, `s=1.609220` vs
  Professional `-230.588` / `1.284286` (target `-230.58 ±0.1`). Before the
  fix `s` stayed at `2.0` and `sigyy=-227.984`. The modification is gated on
  `s >= 1`, so `hypo7`/`hypo9` (no structure record, `hyhis4=0`) are
  numerically unchanged (`.dbs` identical except the nondeterministic
  `time_calculation`). Regression: sfnet 197-test classification identical
  to the a02e7c2 baseline (160 PASS / 37 RUNFAIL / 0 PARSE); hypo1-4 `.dbs`
  byte-identical to both the integrated and the a02e7c2 identity baselines;
  internal suite 16/16 (clean build). Evidence:
  `evidence/structure-eq21/`.
- Corpus vs the Professional (re-measured 2026-09-06 after the hypo_ ABI
  fix of hypoplas.cc; values unchanged by the fix — the deviations are
  genuine kernel-level gaps, not memory contamination): hypo7/8/9
  (masin clay) `sigyy` -224.08/-216.29/-224.08 vs Professional
  -231.81/-230.59/-231.81 (3-6 %, targets tol 0.1, rc=1); hypo12 (masin
  visco JM, probe sin `print_apply -no` para volcar .dbs) `sigyy`
  -282.46 vs Professional -143.49 con "iterative solver broke down" en el
  camino (hyhis0 0.5906 estable); hypo13 (Niemunis visco) aborta con
  "severe error in Masin hypoplasticity" tras desbordar la razón de vacíos
  (hyhis0 ~1.86e14) — calibración visco PENDIENTE (no recalibrada en este
  lote). Los targets Pro de hypo12/13 son -143.495635 (tol 1e-3) y
  -144.49 (tol 1e-2).
- **Diag 2026-09-07 (see hypoplasticity_kernels.md): the hypo7/8/9 corpus
  gap is upstream of the kernel.** Truncation series of the SAME .dat
  (axisym + materi_velocity_integrated + group_materi_memory
  -updated_linear) show the codes accumulate different strain states:
  Professional eptyy = -0.300000000000 exactly (linear strain increments)
  vs GNU eptyy = -0.35657 = -ln(1-0.3) (strain increments evaluated in the
  current geometry -> logarithmic accumulation). The void ratio integrated
  by the masin kernel is consistent with each imposed path (Professional e
  1.03->1.1107, GNU 1.03->1.1263 with a deeper dip 0.941 vs 1.006), and
  the 3-6% sigyy deviation follows from the different paths. The masin.c
  kernel itself is not at fault (validated vs umat_hcea.for to 5e-7 at
  driver level; removing materi_strain_plasti changes nothing). The
  -updated_linear strain accumulation lives upstream of the
  hypoplasticity dispatch and must be aligned there (kinematics owner).

## External dependencies

- `masin.c` is self-contained (only `<math.h>`, `<stdio.h>`, `<stdlib.h>`);
  no f2c, no tochnog headers. The `masin_umat` prototype is declared
  `extern "C"` in `hypoplas.cc`.
- Core `db()`, `db_active_index()`, `pri()`, `array_add`/`array_subtract`
  from tochnog.

## Hardcoded parameters / pending refactorings

- `mtesting=0` is hardcoded; the PLAXIS first-call loose-tolerance mode
  (`testing=1`) and the stiffness-only mode (`testing=2`) are not exposed.
- **Intergranular strain end-to-end discrepancy (P4-B2b)**: the kernel matches
  the Fortran, but the FE result in `hypomasin3.dat` differs (~20%). Root
  cause hypothesis: the equilibrium iterations inside each substep re-integrate
  with the intergranular delta of the previous iteration while keeping the
  stress of the start of the step, damping the small-strain stiffness. A fix
  would pass the converged delta of the substep into the next iteration (or
  store the intergranular tensor in `old_epi`/`new_epi` and use
  `materi_strain_intergranular`). NOT yet done.
- The Niemunis visco law (`Dr Iv`) is implemented from the manual theory
  with clay-derived defaults for `ee0/pe0/lambda/betaR` (the manual exposes
  only `Dr Iv`). There is NO reference UMAT to validate it numerically; the
  regression test pins the behaviour and the creep exponent is clamped for
  stability. If a reference implementation is later found, it should be
  ported and validated the same way as the JM kernel.
- `materi_history_variables >= 8` is enforced with an explicit error; the
  check could be lifted to a softer warning for backwards compatibility.
- The `OCR` initial-void-ratio formula duplicates the Fortran; verify against
  the Fortran when the OCR path is exercised.
- `theta` in `_strain_intergranular_masin_clay` is accepted but unused (the
  kernel interpolates with chi).
- `masin.c` uses fixed `props[29]` and `statev[16]` sizes; refactor to
  structs if more variants (strength reduction, visco) are ported.
- **Basic route**: the structure evolution Eq. (21) is implemented
  (`5a98958`). Intergranular strain is still not available for the basic
  route (the basic kernel has no `istrain` branch). The manual Eq. (12) does
  NOT print the pyknotropy modification of Masin (2007) Eq. (8); it was
  taken from the paper referenced by the manual ([12]) and validated
  end-to-end against `hypo8`.

## Differences with the Professional version

- **hypoplasticity history slots and stress path** (known difference, not
  decided; AUDIT K14). GNU: Leaves `hyhis1..5` = 0 (`hypo7/8/9`); `hypo3`
  -0.1034; `hypo10` disx 4.3e-9. Professional 02-08-2026: Fills the slots
  (e.g. `hyhis3` 22519); `hypo3` -0.1629; `hypo10` 9.84e-8. Status: UNKNOWN
  root cause (partly GNU-suspect); investigate before touching code. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
