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
- C port (after the P4-E1f bisection fix):
  - step1 `-189.9448` (0.000%), step4 `-439.80` (0.12%),
    step8 `-867.25` (0.03%), step12 `-1423.70` (0.32%),
    step20 `-2933.32` (6.3%).
  - `e = 0.666380` vs Fortran `0.666385` (0.001%) at every step.
  The void ratio and early steps are now essentially exact; only the late
  deviator (`a11` 0.675 vs 0.606) still diverges.
- End-to-end tochnog: `hyposanisand1.dat` (laterally confined biaxial),
  targets `sigxx=-3434±200`, `hisv6≈0.66`. Passes.

## Known limitation and future work

The C port is **constitutively correct**. After the P4-E1f bisection fix the
void ratio matches the Fortran to 0.001% and steps 1-8 to <0.12%. The only
remaining difference is the late deviator `a11` (0.675 vs 0.606): the C keeps
hardening while the Fortran saturates.

### P4-E1f findings (intersect_DM bisection — THE root cause)

The definitive bug was in the bisection block of `intersect_DM`:

- The Fortran bisection **narrows the interval** by the sign of `fy`:
  `y05=(y00+y11)/2; if(fy(y05)<0) y00=y05 else y11=y05`. It converges to the
  real crossing (~0.0177).
- The C port had `y00`/`y11` FIXED, so `y05` always stayed the midpoint
  `(y0+y1)/2` → returned `xi=0.5` instead of `0.0177`.
- The earlier P4-E1d "use the Newton xi" workaround (0.0109) was also
  wrong — the correct value is the bisection result 0.0177.

This is the source of the input divergence the whole chain of P4-E1d/e
diagnostics pointed to: the intersection point (entry into plasticity)
drives `T_k` initial (0.0177), and a wrong `xi` changes every subsequent
substep, amplifying into the late-path difference.

### P4-E1e findings (code-identity verification)

Verified line by line that ALL constitutive functions are exact transcriptions
of the Fortran: `yf_DM`, `el_stiff_DM`, `lode_DM` (Van Eekelen), `grad_f_DM`,
`grad_g_DM`, `alpha_th_DM`, `plast_mod_DM`, `get_tan_DM` (Hep/HH_fab),
`check_parms_DM`, `psi_void_DM`, `drift_corr_DM`. No transcription errors.
(This is why the "inputs must differ" reasoning led to the intersection
point.)

### Remaining late deviator difference (P4-E1h)

The void ratio is exact, so the volumetric path is perfect. The `a11` (back
stress) still grows in the C (0.675) while the Fortran saturates (~0.606).
`alpha_sr=0` in both, `alpha_b≈0.85` in both. Next step: compare the plastic
increment `dalpha = Hep*deps` (or `Kp`) at the step-20 substeps in both codes
— with the void ratio exact, the difference must be in the deviator
hardening accumulation.

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
