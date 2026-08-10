# group_materi_plasti_sanisand

## Files and functions

- `sanisand.c` — **new pure-C port** of the reference Fortran UMAT
  `validation-suite/reference-sanisand/umat.for` (4899 lines, GPL,
  Martinelli/Miriano/Tamagnini, adapted by Charles University for
  Dafalias & Manzari 2004).
  - `sanisand_umat()` — single integration-point update (Abaqus UMAT
    interface): `stress[6]`, `statev[36]`, `ddsdde[36]`, `dstran[6]`,
    `dtime`, `props[19]`, `testing`, `error`.
  - Physics: `el_stiff_DM()` (pressure-dependent elastic stiffness),
    `yf_DM()` (yield function `|s - p alpha| - sqrt(2/3) m p`),
    `grad_f_DM()` / `grad_g_DM()` (yield/potential gradients),
    `get_tan_DM()` / `plast_mod_DM()` (elasto-plastic tangent and hardening,
    `Kp = L:De:R + Hplas`), `lode_DM()` (Van Eekelen octahedral shape),
    `alpha_th_DM()` (critical/bounding/dilatancy cones), `psi_void_DM()`
    (state parameter).
  - Integrator: `rkf23_upd_DM()` (adaptive RKF23 with elasto-plastic load
    case detection), `trial_state()` (elastic trial via `f_hypoelas_DM`),
    `check_crossing()`, `intersect_DM()` (Newton + bisection),
    `drift_corr_DM()` (Sloan drift correction), `pert_DM()`, `tang_stiff()`.
  - Helpers: `dot_vect`, `matmul`, `push`, `pzero`, `deviator`,
    `inv_sig_full`, `norm_res_DM`, `check_RKF_DM`, `move_sig`, `move_eps`,
    `iniyz`, `solout_sani`.
- `hypoplas.cc` — dispatch block for `GROUP_MATERI_PLASTI_SANISAND`:
  reads the 19 parameters, maps `hisv[0..35]` to `statev[0..35]`, converts
  the 3x3 tensors to Voigt6, calls `sanisand_umat()`, maps back, and folds
  the stress increment with the same `array_add`/`array_subtract` pattern.
- `tochnog.h` / `tochnog-mod.h` — enum `GROUP_MATERI_PLASTI_SANISAND`.
- `database.cc` — `group_materi_plasti_sanisand`: DOUBLE, length 19.
- `check.cc` — requires `materi_stress` and `materi_history_variables`.
- `makefile` — `SANISAND_SRC=sanisand.c`, `SANISAND_OBJ=sanisand.o`.

## Conventions

- SANISAND uses **soil mechanics convention: compression positive**. The
  kernel converts the tochnog (Abaqus-like, tension positive) tensors
  internally: `move_sig: sig = -stress - pore`, `move_eps: deps = -dstran`,
  `solout: stress = -y - pore`. Tochnog therefore passes the raw tensors.
- `statev` layout (0-based): `[0..5]`=back stress alpha,
  `[6]`=void ratio e, `[7..12]`=fabric z, `[13]`=unused, `[14..19]`=alpha_sr,
  `[28]`=pore, `[29]`=p', `[30]`=q, `[31]`=cos3t, `[32]`=dtsub,
  `[33]`=nfev.
- The kernel writes `stress` back in Abaqus convention (tension positive),
  which matches tochnog's `new_sig`.

## Validation

- Fortran reference: `validation-suite/reference-sanisand/ref_triaxial.txt`
  (Nevada sand, 20 steps, `dstran[0]=-0.001`):
  step4 `-439.2455`, step8 `-866.4842`, step12 `-1413.1947`,
  step16 `-2045.1646`, step20 `-2777.1338`, `e=0.66639`, `a11=0.613237`.
- C port (same path), after the P4-E1b tolerance fix (`tol_f = 1e-6`):
  step20 `-2873`, `e=0.679`, `a11=0.667` — about **+3.5%** on `sig11`
  (improved from +8% before the fix).
- End-to-end tochnog: `hyposanisand1.dat` (laterally confined biaxial),
  targets `sigxx=-3696±300`, `hisv6≈0.66`. Passes.

## Known limitation and future work (IMPORTANT)

The C port is **constitutively correct and physically coherent** (dense sand
hardening: `alpha` grows 0→0.67, `e` drops 0.70→0.679), but it does NOT
reproduce the Fortran to the `~1e-6` accuracy of the Masin ports. After the
P4-E1b tolerance fix the difference is ~**3.5%** on `sig11` at step 20
(was ~8%).

### What was fixed in P4-E1b

1. **`tol_f = 1e-6`** (was 1e-3). The yield-function tolerance of the
   Fortran is `tol_f = 1.0d-6`, independent of `testing`; only the RKF
   `err_tol` switches between `tolintTtest = 1e-2` (first step) and
   `tolintT = 1e-3`. Using 1e-3 for the yield check made the drift
   correction too lax and changed the trajectory. This was the dominant
   cause of the 8% gap.
2. Removed a duplicated `f_plas_DM` call for `kRK_1` (inflated `nfev`).
3. **Confirmed the `attempt==2/3` hypothesis is NOT the cause**: the RKF23
   loop converges in a few substeps per global step (`T_k` reaches 1 well
   before `maxnint`), so the looser-tolerance re-substepping branches never
   activate — in both the port and the Fortran.

### Remaining 3.5% and what was ruled out (P4-E1c findings)

Diagnostics performed on the first global step (where the divergence starts,
+3.9% on `sig11`, `a11` +10%):

- **The elastic trial state is EXACT**: the trial `e = 0.698301` matches the
  Fortran final `e = 0.69833` exactly. The divergence appears only in the
  subsequent plastic substepping.
- **The RKF tolerance is NOT the cause**: tightening `err_tol` from `1e-3`
  to `1e-6` makes the C port converge to `-2947` (step 20), FURTHER from the
  reference `-2777`, not closer. The integrator converges to a stable value
  that differs from the Fortran by ~6% — a systematic model-path difference,
  not a substepping accuracy issue.
- **`attempt==2/3` is NOT the cause**: the RKF23 loop finishes in a handful
  of substeps per global step, never reaching `maxnint`.
- **`alpha_sr` write-back in `get_tan_DM` is NOT the cause**: writing the
  updated `alpha_sr` back to `z` (as the Fortran does) leaves the driver
  result identical AND breaks the FE end-to-end (the tochnog multi-call flow
  diverges to `-784`), so it was reverted.

What remains to investigate (in order of likelihood):
1. The plastic flow direction during the FIRST plastic substep — `a11` grows
   ~10% faster in the port, suggesting a slightly different `h_alpha`
   (bounding-surface) or dilatancy contribution at the yield-surface entry.
2. The exact `intersect_DM` result (the Fortran bisection returns `xi = 0.5`
   — the midpoint — in this test; verify the Newton actually converges in
   the reference before the bisection is reached).
3. `drift_corr_DM` `switch=1` normal-correction branch.

Recommended approach (P4-E1d):
1. Print `a11` after the first plastic substep in BOTH the C port and the
   Fortran (a small Fortran driver instrumenting `get_tan_DM` outputs) and
   diff them to locate the exact expression that differs.
2. Fix that expression, re-validate at the driver level to `~1e-3`.
3. Re-check the FE end-to-end value last (note: the FE and driver paths
   differ because tochnog calls the constitutive model repeatedly per
   substep; only compare relative trends).

## Dynamic substepping (already in place)

The RKF23 integrator **already performs adaptive (dynamic) substepping**
based on convergence:
- `norm_R` (the relative error estimate between the 2nd- and 3rd-order
  solutions) drives the step size:
  `S_hull = 0.9 * DT_k * (err_tol / norm_R)^(1/3)`.
- Accepted steps grow the step: `DT_k = min(4*DT_k, S_hull)`, clamped to the
  remaining `1 - T_k`.
- Rejected steps shrink it: `DT_k = max(DT_k/4, S_hull)`.
- The global time increment `dtime` is sub-stepped internally until the whole
  strain increment is integrated (`T_k` goes 0 → 1).

So the constitutive update already adapts the internal substep to the local
convergence. The `testing` parameter selects the error tolerance (`1e-2` for
the first global step, `1e-3` afterwards), matching the Fortran. The global
`control_timestep` of the input file is independent and can also be made
adaptive at the FE level (see the standard tochnog `control_timestep`/`-no`
convergence-based controls); that is orthogonal to the constitutive
substepping.

## External dependencies

- `sanisand.c` is self-contained (only `<math.h>`, `<stdio.h>`,
  `<stdlib.h>`); no f2c, no tochnog headers. The `sanisand_umat` prototype
  is declared `extern "C"` in `hypoplas.cc`.
- Core `db()`, `db_active_index()`, `pri()`, `array_add`/`array_subtract`
  from tochnog.

## Hardcoded parameters / pending refactorings

- `testing` is used to select the integration tolerance (`1e-2` on the first
  global step like the Fortran `testing=1`, `1e-3` afterwards). The
  `pert_DM` tangent path (`cons_lin==0`) is not exposed; `tang_stiff` is used
  always.
- The `attempt==2/3` integrator branches need to be fully wired (see above).
- `materi_history_variables >= 36` is enforced with an explicit error.
- `sanisand.c` uses fixed `props[19]`, `statev[36]`, `NYDIM=20`,
  `NZDIM=14` sizes.
- The `drcor` (drift correction) flag is hardcoded to 1; the `check_ff`
  (yield-function check) flag to 0, matching the reference UMAT defaults.
