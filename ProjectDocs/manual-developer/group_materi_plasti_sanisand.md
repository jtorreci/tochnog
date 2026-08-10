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
- C port (same path): step20 `-3023`, `e=0.671`, `a11=0.660` — about
  **+8%** on `sig11` (grows from +1% at step 4).
- End-to-end tochnog: `hyposanisand1.dat` (laterally confined biaxial),
  targets `sigxx=-3638±300`, `hisv6≈0.66`. Passes.

## Known limitation and future work (IMPORTANT)

The C port is **constitutively correct and physically coherent** (dense sand
hardening: `alpha` grows 0→0.66, `e` drops 0.70→0.671), but it does NOT
reproduce the Fortran to the `~1e-6` accuracy of the Masin ports. The
difference (~1% at step 4 growing to ~8% at step 20) comes from the
**adaptive elasto-plastic substepping**, not from the constitutive equations.

### Root cause

The reference integrator `rkf23_upd_DM` has a three-level acceptance logic:

1. `attempt == 1`: strict tolerance `tolintT = 1e-3`, `maxnint = 50000`.
2. `attempt == 2`: when `maxnint` is exceeded (or the drift flag `switch3`
   is set), the tolerance is relaxed to `err_tol_1 = 1000 * tolintT` and
   `maxnint_1 = 2 * maxnint`, then the step is retried.
3. `attempt == 3`: final "accept the solution" fallback with drift
   correction on `y_k`.

In this port the `attempt == 2/3` branches are present but the conditions
that activate them (the exact `switch3`/`mario2`/`ksubst > maxnint_1`
triggers) were simplified, so the looser-tolerance re-substepping never
actually runs. As a result the substepping path differs from the reference
and, because the elasto-plastic model is sensitive to the strain-path
subdivision, the stress drifts by up to ~8%.

### Recommended path to close the gap (future work, P4-E1b)

1. **Replicate the full `attempt` state machine**: make `attempt`,
   `maxnint_1`, `err_tol_1`, `err_tol_n`, `mario2` and `switch3` follow the
   Fortran exactly — in particular the block that triggers `attempt=2` when
   `(ksubst > maxnint_1) || (switch3 == 1)`, and the `attempt==3` fallback
   at the bottom of the loop.
2. **Propagate the `plastic` flag by reference** like the Fortran common
   block: it is shared between `get_tan`, `plast_mod`, `drift_corr` and the
   main loop. Currently it is passed by value in several call sites, which
   can desynchronise the elastic/plastic decision.
3. **Re-validate** against `ref_triaxial.txt` after each change, targeting
   `~1e-3` on `sig11` at step 20 (the reference prints 4 decimals).
4. Only after the driver-level match, re-check the FE end-to-end value in
   `hyposanisand1.dat`.

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
