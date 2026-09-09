# DEV-MASIN3: response-curve fitting of the Professional masin kernel

Sprint 2026-09-09, worktree `/tmp/opencode/tn-dev-masin3` (branch
`dev/masin3`, from 40fa51e). Method: systematic single-element probes at
controlled initial states run in BOTH binaries with identical fully
imposed nodal velocities (fresh dir per run), cross-checked with a
standalone kernel harness (`masin.c` compiled outside the FE) and a
1-element trajectory replay with a Newton solve on the lateral strain.

All probe data: `/tmp/masin3_validation/` (`probes/`, `analysis/`,
`scratch/{kprobe,replay}`, `corpus*` runs). Generator:
`gen_probe.py`; runner `run_tn.sh`.

## Headline result

**For `group_materi_plasti_hypo_masin_clay` the Professional's kernel
stress response IS the public hcea formulation** — measured identical
(<= 1e-4 relative, i.e. at the FE print precision) across:

- isotropic and K0-ish initial states (p = 50/78.6/113.8/144.5 kPa),
- imposed lateral increments depsx = 0 .. 2e-5 with depsy = -1e-5,
- plane strain and axisymmetric, dt = 10 (corpus) and 1e-4,
- with and without the intergranular record (istrain = 0 and 1),
- 10x smaller and 10x larger increments (ratio stable to ~1%).

Example (probes `pl0_x*`, plane strain, sig = -50 iso, e = 0.6,
depsy = -1e-5): Pro and GNU dsigxx/dsigyy agree to the printed digit at
every depsx; at depsx = 0 both give exactly (-0.00529, -0.046537).

**The corpus residuals (hypo7/8/9/12/13) are NOT a kernel-level residue
in the clay family.** They are the product of the Professional's
EQUILIBRIUM-LOOP semantics (tangent + iteration behaviour), its
STATE-FEED details, and — for the basic `-masin` record — a different
kernel. Evidence per mechanism:

## 1. The FE-loop path difference (hypo12/13)

GNU runs exactly 2 equilibrium iterations per step (top.cc default
`iteration_min = iteration_max = 2`), assembling the element stiffness
from the kernel's returned `ddsdde` (`perturbate()` = `fs*LL` for
istrain=0, `m_R*LL` for istrain=1). Measured fingerprints:

- **hypo13 (istrain=0)**: the Professional's lateral strain per step is
  FROZEN at exactly `depsx = -0.25 * depsy = +2.5e-6` for 1000+ steps
  (corpus13_p `eptxx`: 2.500e-6 at steps 1..1000, then drifting to
  1.09e-5 late) while GNU converges to the true lateral-neutral point
  (0.96e-6/step; kernel-map neutral = 0.95e-6). The ratio 0.25 equals
  the plane-strain elastic prediction `nu/(1-nu)` with `nu = nu_pp =
  0.2` only if the tangent carries normal-normal Poisson coupling —
  the public hcea `LL` tensor is diagonal in the normal components
  (`LL(i,j)=0` for i!=j normals), which predicts ~0.10 (GNU's value).
  The Professional's loop therefore behaves like a one-pass/loose
  solve with an ELASTIC isotropic tangent (frozen elastic prediction),
  not like the converged Newton GNU performs.
- **hypo12 (istrain=1)**: both binaries evolve (`depsx` 2.4982e-5 Pro
  vs 2.5391e-5 GNU at step 1, both growing similarly): tangent ratio
  Pro/GNU = 0.984 (Pro 1.6% stiffer here). The per-step dsigyy differs
  by ~0.9%, accumulating to the 0.101 kPa final gap.
- Scalar scaling of the returned `ddsdde` (test factor `M3_TAN` in a
  scratch build) reproduces the Pro's step-1 lateral path exactly
  (`M3_TAN=0.38444` gives epsx = 2.5000362e-6 vs Pro 2.4999933e-6) but
  destabilizes the full 10000-step run (epsx -> 1.5e8): a CONSTANT
  scalar is the wrong functional form — the Professional's loop
  semantics cannot be emulated by a kernel-side tangent scalar.

## 2. State-feed difference: the sensitivity slot (hypo12/13)

The Professional feeds `hyhis4` (the sensitivity/structure slot) into
the kernel **even when no structure record exists**; its step-1 `.dbs`
dump shows `hyhis4 = 2.0` (the corpus `control_reset -hyhis4 2`) while
GNU gates the feed to `struct_active` and feeds 0. With `sensit = 2`
and `kparam = 0` the hcea kernel takes the softmodel branch with
`Stf = 1` (no stiffness change) and `N* = N + lambda* ln(2) = 1.0693`,
which doubles `peast`, halves `fd`, and halves the asymptotic (NN)
stress contribution: kernel response at the corpus step-1 input
changes from (-0.00529, -0.04654) to (+0.01169, -0.05215), bracketing
the Professional's corpus step-1 stress (-0.0000000000, -0.05157).
Feeding sensit=2 in GNU ALONE (state-feed experiment, since reverted)
made hypo12 worse (-158.67) because GNU's converged path is not the
Professional's path: the state feed and the loop semantics form one
coupled system and cannot be transported piecewise.

## 3. A Professional-side visco-route quirk (hypo13, first step)

The Professional's corpus visco run starts with a void-ratio increment
10x the strain bookkeeping: step-1 dump shows `hyhis0 = 0.599879999893`
(`de = -1.2e-4 = 10 * (1+e) * tr(deps)` with `tr(deps) = -7.5e-6`,
`dtime = 10`) — consistent with `(1+e)*tr(D)*dtime` applied on top of
the increment-form input. The effect is startup-only (the final
`hyhis0 = 0.603664` is again consistent with the plain bookkeeping of
the whole trajectory). Not reproduced in GNU; documented.

## 4. The basic `-masin` record runs a DIFFERENT kernel in the
   Professional (hypo7/8/9)

Imposed-strain probes at the hypo7 state (sig = -100 iso, e = 1.03,
phic = 22.56, lambda* = 0.11, kappa* = 0.016, N = 1.375, r = 0.4):

| probe (dsigyy, kPa)     | Professional | GNU hcea | ratio |
|--------------------------|--------------|----------|-------|
| oedometric (eps=(6.3e-5,-5e-4,6.3e-5)) | -2.562 | -1.429 | 1.79 |
| 10x smaller increment    | -0.2556      | -0.1426  | 1.79 |
| 10x larger increment     | -26.93       | -14.67   | 1.84 |
| unload (eps_y = +5e-4)   | +3.501       | +2.378   | 1.47 |
| constant-p path          | -1.677       | -0.178   | 9.4  |
| K0 start (sig=-50/-100/-50) | -2.378    | -1.122   | 2.12 |
| p = 25 / 200 (oed)       | 1.63 / 2.27  | -        | -    |

The Professional DOES read `r` (r = 0.2/0.4/0.6 give dsigyy =
-4.248/-2.562/-1.995) but the response is not the hcea kernel with any
single substituted parameter: at ISOTROPIC states every path and
increment is reproduced by the hcea kernel with an effective Poisson
ratio `nu_eff` (r = 0.2 -> nu_eff = 0.054; 0.4 -> 0.243; 0.6 -> 0.315,
fitted from dsigxx AND dsigyy simultaneously), but at the K0 initial
state no `nu` matches (Professional is 60% stiffer in dsigyy at the
K0 state). The volumetric response matches everywhere (r_v = 1.00,
probe table "r_v" column). Conclusion: the Professional's basic
`-masin` is an isotropic, r-parameterised formulation that coincides
with hcea(nu_eff(r)) at isotropic states — consistent with an
older Masin-clay implementation, NOT with parameter substitution into
hcea. Closing hypo7/9/8 needs that formulation, not a constant delta.

## 5. Falsified candidates (this sprint)

- Constant stiffness multiplier: ratio varies 1.63..2.34 with p/e
  (hypo7 family) and is 1.000 for the clay family — refuted.
- RKF tolerance / substepping: step-1 ratio invariant over 100x
  increment sizes (1.792/1.793/1.836) — refuted as the driver.
- alphaG / alphaE / alphanu anisotropy defaults (sweeps 0.25..2.0):
  wrong signature (Pro responses always have dsigxx = dsigzz;
  anisotropic families break that) — refuted.
- r into the G0 slot (intergranular activation): m_R ~ 0 response —
  refuted.
- Strain-measure feeding alone (log vs engineering): a validated
  trajectory replay (`scratch/replay`, reproduces GNU's corpus hypo13
  final to 7e-4 kPa) shows engineering feed alone moves hypo13 to
  -146.05, not -144.50 — refuted as the sole mechanism (it IS a real
  secondary difference: GNU accumulates log strain — printed eptyy
  -0.3567 = ln(1/0.7) — while the Professional stores engineering
  strain -0.300; the kernel therefore receives different per-step
  increments).
- 2020 hcea variant / porting bugs in masin.c: masin.c is
  line-verified against umat_hcea.for (incl. `LL(i,j)=LLfour(i,i,j,j)`)
  and the clay-family probes return identity — refuted.

## Per-test state after this sprint (source unchanged, documented)

| test | GNU (before = after) | Professional | target | gap |
|------|----------------------|--------------|--------|-----|
| hypo7 | -230.792534395 | -231.808976806 | -231.80+-0.1 | 0.44% |
| hypo8 | -217.670847523 | -230.588113239 | -230.58+-0.1 | 5.6% |
| hypo9 | -230.792534395 | -231.810894524 | -231.81+-0.1 | 0.44% |
| hypo12 | -143.393446689 | -143.494849517 | -143.4956+-0.001 | 0.07% |
| hypo13 | -145.690839954 | -144.495200493 | -144.49+-0.01 | 0.83% |

Byte-neutrality of the final build re-verified against the pre-sprint
binaries: hypo7/8/9/12/13 `.his` byte-identical; suite 16/16.

## What closing the targets requires (next sprint)

1. hypo12/13: reproduce the Professional's equilibrium-loop semantics
   (loose/two-pass solve with an elastic isotropic tangent) for the
   masin family — `top.cc`/`so.cc`/`calcul*` (currently outside the
   masin-family edit surfaces) or an explicitly gated per-family loop
   switch; a kernel-side emulation is impossible (mechanism above).
2. hypo7/9/8: implement (or fit) the Professional's basic-masin
   formulation (section 4). The imposed-probe surface in
   `/tmp/masin3_validation/probes` is the fitting data set.
3. Optional fidelity extras: sensit=2 state feed (with the loop
   semantics, not alone), visco first-step void quirk.

# DEV-EQLOOP addendum (2026-09-09): the loop hypothesis falsified, mechanisms relocated

Follow-up sprint (`dev/eqloop`, same base 40fa51e) with solver-layer
edit permission. Full evidence: `eqloop_semantics.md` (same directory);
probe matrix in `/tmp/eqloop_validation/probes/`.

1. **Section 1's "equilibrium-loop semantics" mechanism is
   FALSIFIED.** The Professional's default loop is 2 passes per step
   (proved by `control_timestep_iterations` A/B: N=1/2/3 give ratio
   1.9e-4/1.7e-6/1.4e-8 and finals -143.49497/-143.49485/-143.49482 on
   hypo12 — default == N=2) with the kernel-returned consistent tangent;
   GNU's 2-pass loop is already at its fixed point (N=2 == N=20 on both
   corpus tests). A converged fixed point is tangent-independent: no
   loop-side change can move GNU toward the Professional.
2. **Section 2's state-feed mechanism is FALSIFIED.** The Professional's
   stress response on imposed paths is IDENTICAL (10 digits) with
   `hyhis4` reset to 0 vs 2 — the stored sensit value does not reach the
   kernel response without a structure record. Item 3 of "What closing
   the targets requires" is void.
3. **The real mechanisms are kernel-route differences**, measured with
   fully imposed strain paths (no free dofs): (a) the `_visco` route is
   ACTIVE in the Professional at corpus parameters — an exactly
   isotropic stress rate -0.008732 kPa/step at the corpus state
   (dtime-exact, Dr/Iv-independent at 1e-6..2e-6, ~p²/peast scaling),
   which pins dsigxx=0 on the elastic-poisson lateral path and produces
   the section-1 "frozen prediction" fingerprint; GNU's visco route is
   inert. (b) the intergranular `istrain=1` route differs by ~1%/step in
   dsigyy at the corpus12 state in BOTH geometries.
4. **DEV-MASIN2 correction:** its "corpus-param Pro visco run
   bit-identical to no-visco" entry is contradicted by the archived
   `corpus13_p/run.log` itself (step-1 dsigyy -0.0515678 = the
   active-visco value, not the no-visco -0.0428).
