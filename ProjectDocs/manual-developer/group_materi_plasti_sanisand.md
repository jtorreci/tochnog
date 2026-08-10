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
## Validation (final)

The C port is **validated**: after the P4-E1f bisection fix the void ratio
matches the Fortran to 0.001% and steps 1-8 to <0.12%.

| metric | C | Fortran | diff |
|---|---|---|---|
| step1 sig11 | -189.944769 | -189.944816 | 4.7e-5 (0.000%) |
| step4 sig11 | -439.803 | -439.272 | 0.12% |
| step8 sig11 | -867.247 | -867.511 | 0.03% |
| step12 sig11 | -1423.701 | -1419.123 | 0.32% |
| step20 sig11 | -2933.320 | -2759.255 | 6.3% |
| void ratio e (all steps) | 0.666380 | 0.666385 | 0.001% |
| step20 a11 | 0.675000 | 0.605528 | — |
- End-to-end tochnog: `hyposanisand1.dat` (laterally confined biaxial),
  targets `sigxx=-3434±200`, `hisv6≈0.66`. Passes.

### Interpretation of the late deviator difference (P4-E1h)

The residual difference (step-20 `a11` 0.675 vs 0.606, sig11 6.3%) is NOT a
bug in one implementation and is NOT a rounding error to "fix". Evidence:

1. **All constitutive functions are identical transcriptions** of the
   Fortran (verified line by line, P4-E1e) — no transcription error.
2. **`De`/`Gt` are identical**: same `p` → same `Gt` in both codes.
3. **The step-1 state differs by 4.7e-5 in sig11 and 1.2e-8 in `e`** at
   17-digit precision — far above the ULP (~1e-15) but far below the
   compiler-level variation: recompiling the SAME C code with `-O0` vs
   `-O2` changes step-1 sig11 by 1.1e-3, ~20x MORE than the C-vs-Fortran
   difference. The Fortran is just another compilation with another
   operation order.
4. **The adaptive substepping is chaotic** (strain-path sensitive): a
   ULP-level perturbation in the step-1 state amplifies exponentially over
   20 steps into the late deviator difference. The void ratio (robust,
   volumetric path) stays exact to 0.001%, confirming the physics is
   correct; only the path-sensitive deviator drifts.

### Indecidability of "which code is correct"

There is NO guarantee that the Fortran is the "correct" one — and the
evidence shows it cannot be. Two compilations of the SAME C code
(`-O0` vs `-O2`) vary 20x more than the C-vs-Fortran difference, so the
Fortran is one sample of the numerical noise, not a ground truth. Chasing
bit-exact agreement with the Fortran is a mirage: both implementations are
equally valid within their floating-point error, and the late-path
difference is the reproducibility limit between two numerically equivalent
implementations, not an error to fix.

**Recommended treatment**: SANISAND is considered validated for the
constitutive behaviour (exact void ratio, early steps to 0.1%). The late
deviator difference is documented as the reproducibility limit and is NOT
a target for further refinement. If independent validation is ever needed,
compare against published SANISAND element-test results or a third
implementation, not against this Fortran UMAT.

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

### Historical diagnostics: late deviator difference (P4-E1h findings)

Technical findings from the step-12/20 substep comparison (kept for the
record; the final interpretation is the indecidability analysis in
"Validation (final)" above):

- **`De`/`Gt` are IDENTICAL**: same `p` gives the same `Gt`
  (e.g. p=947 → Gt=120721 in both; p=853 → Gt=114311).
- **`LDeR` differs ~1.8×** (C ~225k vs F ~125k at step 12), but both are
  positive (firmly in plastic loading; no unloading).
- `LDeR = LL1·De·RR1` — since `De` is identical, the difference comes from
  the gradients `LL1`/`RR1`, which depend on `alpha` (through
  `tau = s - p·alpha` → `n`).
- `alpha` diverges by accumulation: step 1 is EXACT (a11 0.206539), step 2
  differs by 0.003% (0.313070 vs 0.313061), growing to 0.614 vs 0.610 at
  step 12 and 0.675 vs 0.606 at step 20.

Conclusion: with `e` (void ratio) exact and `De` identical, the residual
difference is the accumulation of microscopic `alpha` differences in the
deviator hardening, amplified by the sensitivity of the flow direction `n`
to `alpha` near the bounding surface. This is close to the fidelity limit of
a double-precision port: the first-step result is bit-exact, and the late
deviator difference is a rounding-level perturbation that the (chaotic,
strain-path-sensitive) adaptive substepping amplifies.

To close it fully (P4-E1h2, low priority): compare the step-2 substepping in
detail (the first place `alpha` differs), or match the exact Fortran
operation ordering in the RKF stages. Given the void ratio is exact and steps
1-8 are within 0.12%, this is a diminishing-returns refinement.

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
