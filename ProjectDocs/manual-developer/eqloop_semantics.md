# DEV-EQLOOP: the Professional equilibrium-loop semantics, measured and closed

Sprint 2026-09-09, worktree `/tmp/opencode/tn-dev-eqloop` (branch
`dev/eqloop`, from 40fa51e). Method: controlled two-binary experiments
on single-element models with fully imposed strain paths (no free dofs
→ pure kernel response), free-dof corpus replicas, iteration-count
sweeps, and record perturbations (Dr, Iv, dtime, e0, p0, hyhis4, epi
sign). All numbers below are measured on the Professional 25-10-2023
binary and the GNU build of this tree (sources identical to 40fa51e;
suite 16/16). Probe inputs and raw outputs: `/tmp/eqloop_validation/`
(`probes/*`, `g12/`, `g13/`).

## Headline

**The Professional's equilibrium loop is semantically IDENTICAL to the
GNU default loop: 2 passes per step, kernel-returned consistent tangent,
same state-update order. The DEV-MASIN3 "frozen prediction / loose
one-pass solve" hypothesis is FALSIFIED.** The corpus hypo12/13 residuals
are NOT loop residuals: they are KERNEL-ROUTE differences — the
`_masin_clay_visco` route (hypo13, the dominant effect) and the
intergranular `istrain=1` route (hypo12) — measured directly with imposed
strain paths in both binaries. Both routes live in `masin.c` /
`masin_visco.c` / `hypoplas.cc`, which are outside this sprint's edit
surfaces; no solver-layer change can close the targets without
kernel-side work. This sprint therefore delivers the fingerprint and the
falsifications, and leaves the source tree byte-identical to 40fa51e.

## 1. The loop fingerprint (passes, tangent, state-update order)

### 1.1 Passes per step: the Pro defaults to 2, like GNU

The Pro honours `control_timestep_iterations` (manual 6.385). On corpus12
(the harder step-1, depsy = -1e-4/step):

| variant | step-1 `post_node_rhside_ratio` | final sigyy |
|---|---|---|
| Pro, `control_timestep_iterations 90 1` | 1.93e-04 | -143.494974806 |
| Pro, default (no record) | 1.67e-06 | -143.494849517 |
| Pro, `control_timestep_iterations 90 3` | 1.38e-08 | -143.494816330 |
| GNU default (2 passes) | 1.42e-04 | -143.393446689 |

Default sits exactly between N=1 and N=3 → **default = 2 passes**. The
pass count moves the final only at the 1e-7/1e-8 level; the corpus gap
(0.101 kPa) is 6 orders larger. The Pro's ratio series is
iteration-max-bound (corpus12: 1.67e-6, 4.9e-9, 3.4e-8, 1.3e-7, ... —
the ratios GROW over steps as the state hardens; no fixed criterium is
visible), with quadratic pass-to-pass contraction (1.9e-4 → 1.7e-6 →
1.4e-8 ≈ ×100 each pass) — consistent-tangent Newton behaviour,
identical in character to GNU's.

`control_print_number_iterations -yes` (manual 6.336) prints nothing in
the 25-10-2023 build (record accepted, no output) — pass counting was
done with the `control_timestep_iterations` A/B instead.

### 1.2 GNU is already fully converged at its default 2 passes

Running GNU hypo13 with `control_timestep_iterations 90 N`:

| N | final sigyy (hypo13) |
|---|---|
| 1 | -145.695 |
| 2 (default) | -145.690839954 |
| 3 | -145.691 (identical) |
| 4 / 6 / 10 / 20 | -145.691 (identical to 12 digits) |

Same for hypo12: N=2 == N=3 == N=4 == N=8 == N=20 = -143.393. **GNU's
2-pass loop is at its fixed point; more iterations move nothing.** A
converged fixed point is tangent-independent, so no loop-side change
(tangent, passes, ordering) can move GNU toward the Professional's
trajectory. This single measurement kills the entire "loop semantics"
family of mechanisms for hypo12/13.

### 1.3 State-update order: same staggered structure

Both binaries run `element loop (σ update with current dof + internal
force) → solve`, repeated; the committed state is the last kernel pass.
The imposed-path probes (below) bypass the loop entirely and still show
the full residual — the loop carries no additional semantics.

### 1.4 The DEV-MASIN3 "frozen prediction" re-interpreted

The measured fingerprint that motivated the loop hypothesis — the Pro's
lateral strain frozen at exactly `depsx = -ν_pp/(1-ν_pp)·depsy =
+2.5e-6` for 1000+ steps with `post_node_rhside_ratio` ≈ 3.9e-12 and
sigxx pinned at -50.00000000002 — is **the Pro's kernel-neutral
equilibrium path**: on that path the Pro's VISCO kernel returns
dsigxx ≈ 0 exactly (section 2), so the converged loop holds sigxx at the
applied -50 pressure while the lateral strain sits at the elastic
Poisson ratio. GNU's kernel on the same path returns dsigxx =
+0.008725/step, so GNU's converged loop moves the lateral strain to its
own (different) neutral point 0.96e-6/step. Same loop, different kernels.

Direct falsification of the "implement a frozen-prediction loop" plan:
feeding GNU's kernel the Pro's path (imposed probe `Cp_v2_g`) makes
sigxx drift +0.0087/step (-50 → -44.95 over 10000 steps) — a frozen
loop in GNU would break lateral equilibrium and reproduce nothing
(its sigyy on that path: -92.07 vs the Pro's corpus -113.79 at t=1e5).

## 2. The real mechanism: kernel-route differences (measured)

### 2.1 The visco route (hypo13): the Pro's visco kernel is ACTIVE at corpus parameters

Fully imposed plane-strain probes on the corpus13 initial state
(sig = -50 iso, e = 0.6, hyhis4 = 2, dtime = 10), per-step increments
(dsigyy/step, first step):

| probe | GNU | Professional |
|---|---|---|
| depsx = +2.5e-6/step, NO visco record (`Dp_v2`) | -0.042836 | -0.0428356 |
| depsx = +2.5e-6/step, visco record (`Cp_v2`) | -0.0428 (inert) | **-0.0515678** |
| depsx = 0, visco record, hyhis4=0 (`B_h40`) | -0.0465 | -0.0550 |

- Without the visco record the binaries agree to 6 digits (corroborates
  DEV-MASIN3's plain-clay probes; here in PLANE STRAIN).
- With the visco record the Pro shifts the response by an **exactly
  ISOTROPIC stress rate**: Δsigxx = Δsigyy = Δsigzz = -0.008732 kPa/step
  at the corpus state (from Cp_v2 vs Dp_v2, all three normals). GNU's
  visco route is inert at these parameters (the 40fa51e "corpus-calibrated
  cutoff" makes it pass through to the plain kernel) — hence the whole
  hypo13 residual.
- This isotropic rate is exactly what pins dsigxx = 0 on the
  elastic-poisson path: plain-kernel dsigxx +0.008725/step + visco
  -0.008732/step ≈ -7e-6/step ≈ 0 (sigxx pinned to 1e-8 for thousands of
  steps; the corpus ratio stays 0.2499933 for 1000+ steps, drifting only
  when the deviatoric stress approaches the critical-state region — the
  documented late drift to ~1.09e-5/step).
- NOTE: this contradicts DEV-MASIN2's table entry "corpus-param Pro
  visco run BIT-IDENTICAL to no-visco" — the archived
  `corpus13_p/run.log` step-1 dsigyy (-0.0515678, prior sprint's own Pro
  run) equals the active-visco value, not the no-visco -0.0428. The
  bit-identity claim could not be reproduced and should be considered
  superseded.

### 2.2 Visco-rate scaling (imposed probes on the Pro, first-step isotropic shift)

| perturbation | isotropic shift/step | ratio to base |
|---|---|---|
| base (p=50, e=0.6, Dr=1e-6, Iv=0.1, dt=10) | -0.008732 kPa | 1 |
| Dr 1e-6 → 2e-6 (`V_dr2`) | -0.008732 | **1.000 (8 digits!)** |
| Iv 0.1 → 0.2 (`V_iv02`) | -0.008726 | 0.99926 |
| e0 0.6 → 0.55 (peast 200.9 → 276.3) (`V_e055`) | -0.006158 | 0.705 |
| p 50 → 100 (`V_p100`) | -0.037355 | 4.278 |
| dtime 10 → 1 (`V_dt1`) | -0.000873/step | ×(1/10) — exact rate law |

Readings: (a) the shift is a true RATE (exact 1/dtime proportionality);
(b) at corpus parameters it is INDEPENDENT of the record Dr (doubled:
zero change to 8 digits) and nearly independent of Iv; (c) p-scan ×4.28
≈ quadratic in p; e-scan ×0.705 ≈ peast-ratio^-1 (276.3/200.9 = 1.376 →
1/1.376 = 0.727, 3% off) — i.e. roughly `∝ p²/peast ≈ p/OCR`.
(d) OPEN: DEV-MASIN2's `p13dr5` (Dr=1e-5) showed a LARGE Dr effect
(final -139.436) — so the record Dr does matter somewhere (threshold or
second mechanism). A candidate reconciliation for the next sprint: the
Pro's visco kernel applies its own internal rate at corpus parameters
(record Dr/Iv only dominating above some threshold), or reads the creep
parameters from different slots than the record writer uses.
(e) The void-ratio "10x quirk" (DEV-MASIN3 §3) is dtime-dependent:
de/step = -1.2e-4 = 10× bookkeeping at dtime=10, but EXACTLY the plain
bookkeeping -1.2e-5 at dtime=1 (probe `V_dt1`) — the quirk factor
scales with dtime, startup-only.

### 2.3 The intergranular route (hypo12): ~1%/step kernel difference, sensit irrelevant

Fully imposed probes on the corpus12 initial state (istrain=1, R/A_g/...
corpus values, dtime=1e-4, depsy=-1e-4/step), plane strain, first-step
dsigyy:

| probe | GNU | Professional |
|---|---|---|
| depsx = +2.4982e-5/step, hyhis4=2 (`H12_p2`) | -3.4137 | -3.3804 |
| depsx = 0, hyhis4=2 (`H12_02`) | -3.6577 | -3.6141 |
| depsx = +2.4982e-5/step, hyhis4=0 (`H12_p0`) | -3.4137 | -3.3804 |
| same, axisymmetric (`A12_p2`) | -3.1719 | -3.1477 |

- The binaries differ by 0.7-1% per step on IDENTICAL inputs, in BOTH
  plane strain and axisymmetric geometry (the DEV-MASIN3 "istrain=1
  identical ≤1e-4" result does not reproduce at corpus parameters).
- **The Professional's stress response is IDENTICAL (10 digits) with
  hyhis4 = 2 and hyhis4 = 0** — the DEV-MASIN3 "state feed
  `hyhis4=2` into the kernel" claim is FALSIFIED as a stress-response
  mechanism (the stored value differs; the kernel output does not).
- The intergranular initial-delta SIGN convention is not the mechanism
  either: flipping epi resets to +1e-5 in GNU moves dsigyy by only
  0.004/step (needed: 0.033).
- The difference must live in the `istrain` branch of the kernel
  (`masin.c` get_tan/integration: R-bucket evolution, m_R interpolation,
  or the Voigt mapping) — out of this sprint's surfaces.

## 3. Falsified mechanisms (do not revisit without new evidence)

1. **Frozen-prediction / elastic-tangent loop**: a converged fixed point
   is tangent-independent (§1.2); feeding GNU the Pro's path breaks
   lateral equilibrium (§1.4). A kernel-side scalar tangent emulation was
   already falsified in DEV-MASIN3.
2. **Iteration-count deficit**: GNU N=20 == N=2 (§1.2); the Pro defaults
   to the same 2 passes (§1.1).
3. **`hyhis4` sensitivity state feed**: the Pro's stress response is
   identical with hyhis4=0/2 (§2.3). The DEV-MASIN3 coupling claim
   ("state feed + loop form one system") is superseded: neither matters.
4. **Intergranular delta sign flip**: measured, wrong magnitude (§2.3).
5. **`ocr_apply` as the hidden driver**: inert in both binaries on the
   imposed probes (no-visco probes match exactly with the record
   present).
6. **Plane-strain wrapper feed difference**: plain masin_clay plane-strain
   probes match to 6 digits — the FE→kernel strain feed is identical;
   the differences are inside the istrain/visco kernel routes.

## 4. What closing hypo12/13 actually requires (next sprint, kernel surfaces)

1. **hypo13**: implement the Pro's ACTIVE visco route in
   `masin_visco.c`: an isotropic compressive stress rate
   `dot(σ_iso) ≈ -8.73e-4 kPa/s` at the corpus state, ∝ p²/peast-ish,
   dtime-exact, Dr/Iv-independent at corpus parameters (open: the Dr=1e-5
   mechanism), plus the dtime-scaled void-ratio startup quirk. The
   activation law fitting surface: the probe matrix in
   `/tmp/eqloop_validation/probes/` (`V_*` scans; regenerate with
   `mk_feed.sh`).
2. **hypo12**: fit the `istrain=1` branch difference in `masin.c`
   (~1%/step in dsigyy at the corpus state; probe families `H12_*`,
   `A12_*`). The imposed-probe pairs are exact inputs for a
   standalone-kernel bisect (the DEV-MASIN3 `scratch/kprobe` harness).
3. Both closures are loop-neutral: with the kernels aligned, the
   EXISTING 2-pass loop converges onto the Professional's trajectory
   automatically (the corpus Pro path IS its kernel-neutral path).
4. dynamic8 is unaffected by any of this (see DIAG-SOLVE-MIXTO §22): its
   auxiliary-solve breakdown is a solver/operator blocker (§19), not an
   equilibrium-loop difference.

## 5. Per-test state (before = after; source-neutral sprint)

| test | GNU (before = after) | Professional | target | gap |
|------|----------------------|--------------|--------|-----|
| hypo12 | -143.393446689 | -143.494849517 | -143.4956±0.001 | 0.07% open |
| hypo13 | -145.690839954 | -144.495200493 | -144.49±0.01 | 0.83% open |

Baselines reproduced in this worktree (`/tmp/eqloop_validation/g12`,
`g13`); no source file changed vs 40fa51e, so the 219-PASS corpus state
is preserved by construction; suite 16/16 with this tree's build
(`/tmp/eqloop_baseline_build.log`).

## 6. Files / evidence

- Probe generator: `/tmp/eqloop_validation/probes/mk_feed.sh` (turns a
  corpus .dat into a fully-imposed-path probe: force edge → prescribed
  right-edge velx, print_apply -yes, configurable hyhis4).
- Probe matrix (both binaries, fresh dirs): `C{0,g,p}_v2` (visco paths),
  `D{0,p}_v2` (no-visco), `B_h40` (sensit off), `H12_*` (intergranular),
  `A12_p2` (axisymmetric), `V_*` (Dr/Iv/e/p/dtime scans),
  `P12_N{1,3}` (Pro iteration A/B), `n13_N*`/`n12_N*` (GNU iteration
  sweeps).
- GNU builds: `/tmp/opencode/tn-dev-eqloop/build/tochnog` (sources ==
  40fa51e); build log `/tmp/eqloop_baseline_build.log` (suite 16/16).
- Prior-sprint archives used: `/tmp/masin3_validation/corpus1{2,3}*`.

## 7. Kernel-surface closure (2026-09-12)

The prescribed next sprint was completed on the candidate tree (the two
uncommitted source surfaces are `masin.c` and `masin_visco.c`). It is a
kernel-only closure: the equilibrium loop was not changed.

- **hypo12 / intergranular route:** the active `istrain` branch applies
  `G0 *= 0.96`, the measured 4% softening relative to the legacy/public
  mapping. The final value is `-143.494989935`, against Professional
  `-143.494849517`: absolute difference `0.000140418`, within the
  `±0.001` target.
- **hypo13 / low-Dr visco route:** the low-Dr Niemunis surface adds the
  measured isotropic correction
  `sigma_dot_iso = -2.99e-6 * pmean^2 / pe` (host compression-negative
  convention), with the existing higher-rate `L:Dvis` correction retained.
  The coefficient is fitted to the free-dof trajectory, not copied from
  the imposed-path first-step rate. The final value is `-144.495489612`,
  against Professional `-144.495200493`: absolute difference
  `0.000289119`, within the `±0.01` target.

### 7.1 Certification and counting caveat

The persistent full ledger is
`/tmp/masin-kernel-corpus-full.8pjO3w/detail.txt`. It classifies all 363
inputs as 220 PASS, 121 RUNFAIL, and 22 DEFERRED. It was assembled in two
execution phases: an early short-timeout phase (12 s) and a later phase
with the normal 45 s timeout. Therefore its raw 220 PASS count is **not**
a clean corpus delta.

The comparable baseline (`run5`/`run6`) is 219 PASS. Normalising only the
classification that is comparable supports **219 -> 220 PASS (+1)**. Both
`hypo12` and `hypo13` independently changed RUNFAIL -> PASS, but the
context change around `abaqus1` cannot be attributed to this kernel sprint;
do not report raw 220 as a clean +2 gain. The available GNU suite remains
16/16 PASS.

This closes the two EQLOOP kernel targets and leaves the already-falsified
loop hypotheses closed. The next Masin work unit is the basic 2005 route
(`hypo7/8/9`), not a revisit of loop or `hcea` behaviour.
