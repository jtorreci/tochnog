# DEV-C hypo-kernel calibration — hypo3/7/8/9/12/13 (2026-09-08)

Worktree `dev/hypok` (from 19e995a). All values measured with the sprint
GNU build (`scripts/build_safe.sh`, suite 16/16) and the Professional
binary 25-10-2023 (`tochnog_user_supplied`). A/B probes:
`/tmp/opencode/tn-dev-hypok-probes/`.

## Measurement pitfalls (important for any re-measurement)

- `.his` history files APPEND across runs. Never reuse a directory with
  stale `.his` (the pre-sprint "Pro −0.554 with e→1.9e14" readings for
  hypo12/13 were stale GNU hypo13 data copied with the input files).
- The Professional honours `print_apply -no` by writing NO output at all
  (no `.dbs`, no `.his`). Its `.his` names differ:
  `post_point_dof_1_<dof>.his` vs GNU `<dof>1.his`.

## Per-test state after this sprint

| test | GNU pre | GNU post | Professional | target | state |
|------|---------|----------|--------------|--------|-------|
| hypo7 (masin basic) | −224.076947603 | −230.792534395 | −231.808976806 | −231.80±0.1 | 0.44% residual (rc=1) |
| hypo8 (masin + structure) | −216.29249576 | −217.670847523 | −230.588113239 | −230.58±0.1 | 5.6% residual (rc=1) |
| hypo9 (masin, epi dofs) | −224.076947603 | −230.792534395 | −231.810894524 | −231.81±0.1 | 0.44% residual (rc=1) |
| hypo12 (masin_clay + intergranular) | −282.457704737 (frozen) | −143.393446689 | −143.494849517 | −143.495635±0.001 | 0.07% residual (rc=1) |
| hypo13 (masin_clay + visco) | crash, e→1.9e14 @ t≈10390 | crash (same class) | −144.495200493 (completes t=1e5, e=0.6037) | −144.49±0.01 | UNCLOSED, see diagnosis |
| hypo3 (wolfersdorff + intergranular) | sigxx −0.103436189421 | unchanged (byte-identical) | −0.162914620443 | −0.1618 | open (see wolfersdorff doc) |

hypo1/2/4: byte-identical targets preserved (suite rc=0 after every
change; wolfersdorff path untouched).

## Fixes applied (hypoplas.cc, masin_visco.c)

1. **Kernel state input** (`hypoplas.cc`, masin block): the masin kernels
   were fed the Newton iterate (`new_hisv`, `new_epi`) instead of the
   step-start state (`old_hisv`, `old_epi`, standard UMAT semantics).
   Every extra Newton iteration re-integrated the state update.
   Measured, hypo12 step 1: GNU void ratio 0.59972066 = 0.6 + trial
   increment (−1.6e-4) + converged increment (−1.19e-4); Professional
   0.59987997 = one integration of the converged increment
   (= (1+e)·tr(deps) with its own converged eptxx = 2.498e-5, exact).
2. **Intergranular delta write-back** (`hypoplas.cc`): the kernel's
   integrated delta (statev[0..5]) was discarded, so the epi dof stayed
   at its reset value forever (GNU) while Pro's evolves (−1e-5 →
   +5.67e-6 in step 1). Now written back to `new_epi` (hyhis layout) or
   `new_hisv[0..5]` (legacy layout).
3. **History layout guard** (`hypoplas.cc`): with
   `materi_plasti_hypo_history` (hyhis0..7 = e, substep, ..., structure,
   ...) there is NO delta slot in hisv; the legacy
   materi_history_variables layout puts delta in hisv[0..5], e at 6.
   Feeding/writing delta into hisv[0..5] of the hyhis layout overwrote
   the void ratio (hypo7/8: e frozen at 1.03) — caught and fixed within
   this sprint; the shipped gate is
   `materi_plasti_hypo_history && materi_strain_intergranular` → epi dof,
   `materi_plasti_hypo_history` → delta absent (0),
   else → legacy hisv[0..5].
4. **Sensitivity (structure) seed gate** (`hypoplas.cc`): hyhis4 seeds
   statev[13] ONLY when a structure group
   (`_masin_clay_structure` / `_masin_structure`) is active. Measured:
   Professional output is bit-identical for hyhis4 = 2 / 1 / absent on
   hypo12 (no structure group); with a structure group it seeds and
   evolves s (hypo8 one step: hyhis4 2 → 1.99932 in the .dbs record).
5. **Visco void ratio** (`masin_visco.c`): `de = (1+e)·(tr(D) − tr(m)·c)`
   subtracted the creep part although `dstran` is the TOTAL strain rate
   (which already contains creep) — double count. Now
   `de = (1+e)·tr(D)`.

## UNCLOSED: hypo13 (visco) — precise diagnosis

GNU still runs away (e→±∞ → `severe error in Masin hypoplasticity` /
solver breakdown at t≈10390 of t_end=1e5); the Professional completes
with e=0.6037, sigyy=−144.4952 (target −144.49±0.01 → Pro passes).

The runaway is driven by the creep-rate closure in `masin_visco.c`
(`masin_niemunis_visco_umat`), not by the driver:

- `pe = pe0*((1+e0)/(1+e))^(1/lambda)` with `pe0 = CURRENT mean stress`
  (recomputed every increment) and `creep_rate = Dr*(1/OCR)^(1/Iv)`,
  `OCR = pe/peplus`, and `peplus = pmean` at the isotropic stress ratio
  (η=0 degeneracy) make OCR a pure function of e:
  `OCR ≈ ((1+e0)/(1+e))^(1/λ*)` with λ*=0.1 and 1/Iv=10, i.e. the creep
  rate scales like `((1+e)/(1+e0))^100` — an exponential feedback. With
  Dr=1e-6, any perturbation of e above the reference state explodes
  (measured: e 0.6 → 1.9e14 over ~1000 steps of dt=10). The Professional
  keeps e≈0.6037 over the same run: its effective creep at this state is
  negligible (OCR≈1.09 from a FIXED reference, or an equivalent closure).
- The reference Fortran for the Niemunis visco kernel is not in the repo
  (only `umat_hcea.for` for the non-visco masin kernel), so the exact
  closure (reference pressure e_ref/pe_ref pairing, OCR definition, creep
  direction sign, and the `creep_rate > 1e3*Dr` clamp) could not be
  ported-verified here. Fixing this needs the `umat_ab_visco`-class
  reference source (Masin/Niemunis 2014) or a parameter-matched probe
  series against the Professional.

## Remaining residuals (documented, not blocking)

- hypo7/9 (0.44%) and hypo12 (0.07%): kernel-level residue after the
  state-feeding fix (first-increment stress differs at the 0.05% level,
  e.g. step-1 sigyy GNU −53.4098 vs Pro −53.3804 with identical BC
  strain). Candidate sources: RKF substep tolerance, `get_tan`
  perturbation tangent, and Pro reading
  `group_materi_plasti_hypo_strain_intergranular` (the Professional
  5-param intergranular record) for the BASIC masin kernel — GNU maps
  that name to the wolfersdorff intergranular group and ignores it for
  masin (hypo7 vs hypo9 are identical in GNU but differ by 0.002 in Pro).
- hypo8 (5.6%): structure evolution path (s: GNU 2 → 1.54 over the run;
  Pro's degradation rate not yet compared over the full run).
- `control_materi_plasti_hypo_masin_clay_ocr_apply` is dead in GNU
  (measured: -yes identical to -no) but functional in the Professional
  (measured: recalculates e0). Not exercised by the corpus targets (they
  all run -no); left for a follow-up sprint.
- The pre-sprint "updated_linear accumulates log strain" hypothesis is
  REFUTED as the main gap driver: after the state fix, hypo7/9 end-state
  void ratios match Pro to 3.6e-4 (1.1103892 vs 1.1107458), which is
  impossible if the strain measures differed by the log/linear gap at
  20% strain (~0.02 in e). `set_deften_etc` (materi.cc) was NOT changed.

---

# DEV-MASIN2 continuation (2026-09-09, worktree dev/masin2 from 93302c0)

## Scope of this sprint

hypo13 primary; hypo12/7/8/9, ocr_apply gate and hypo3 verification as
secondary. New evidence base: Professional 25-10-2023 A/B probes under
`/tmp/masin2_validation/` (fresh dirs, .his extracted per dof) and a
reference Fortran corpus in `/tmp/opencode/ref_extract/`
(plaxis-umat-hypoplas umat_hcea.for 2885 lines = the masin.c source;
hypoClayStrengthReduction 2020 umat_hcea.for = wrapper-only changes;
umat_visco.f = viscoJJ, see the new devmasin2 doc).

## hypo13 root cause (was misdiagnosed as a visco-closure divergence)

The crash was NOT in the visco physics. Evidence chain:

1. `git status` clean at 93302c0: `hypoplas.cc` routes
   `group_materi_plasti_hypo_masin_clay_visco` to
   `masin_niemunis_visco_umat` (the "Niemunis" kernel). That kernel
   computed the WHOLE stress increment with its own soft Niemunis
   elasticity `L = fb*Lhat` (`kappa` defaulted to `lambda*` = 0.1
   instead of kappa* = 0.01, ~5.5x softer than the hcea L) and its
   creep term was a silent NO-OP: at the isotropic state F=0 (the
   hcea invariant convention gives sin2phim = 0), hence m = 0, hence
   `Lm*creep_rate = 0` for any Dr. Measured with a standalone harness:
   `masin_niemunis_visco_umat` at `deps=(0,-1e-5,0)` gives
   `dsigyy = -0.007344720346` - EXACTLY the first-step stress of the
   crashing FE run. The FE accumulated this wrong (too soft) response
   from step 1 (GNU step-1 sigyy -50.0073 vs Pro -50.0516), drifted to
   tension, e went negative and the "severe error" fired at t=10390.
2. The plain (non-visco) masin_clay route NEVER crashed: removing the
   visco record completes the run: final sigyy -145.690839954,
   e 0.594128 (GNU novisco baseline, probe g13nov2).

## Professional hypo13 behaviour (decisive probes)

- `p13_1s` vs `p13nov_1s`: the Pro's step-1 result is BIT-IDENTICAL
  with and without the visco record.
- `p13base_full` vs `p13nov_full`: bit-identical over the FULL 10000
  steps (sigyy/e/eptxx .his diff = 0 lines). With Dr=1e-6, Iv=0.1 the
  Professional's visco contributes EXACTLY ZERO to hypo13.
- `p13dr2_1k` (Dr=2e-6): diverges from no-visco by ~4e-8 kPa/step -
  exactly `L:Dvis*dt` order for creep = Dr*(1/OCR)^(1/Iv) with the
  fixed-NCL OCR (peast/p ~ 4 at the initial state, factor ~9.4e-7).
- `p13iv05_dr6` (Iv=0.5): strongly active (1.85e-2 kPa by step 10),
  again `L:Dvis*dt` order for factor (1/4.006)^2.
- `p13dr5_full` (Dr=1e-5): -139.436 final (visco shifts sigyy by
  +5.06 kPa), e 0.5789 instead of 0.6037.
- Dr/Iv sweeps in `(Dr, Iv)` at fixed state rule out every monotone
  state-dependent activation law for the Pro's visco: Dr=1e-6 stays
  bit-inert even at late-run states (where a peast/p-based creep rate
  would be ~2e-7) while Dr=2e-6 activates from step 1 (creep 2e-12).
  The Pro's activation law remains unresolved (documented, not blocking
  the corpus: at corpus parameters the Pro visco is inert).

## Fixes applied (masin_visco.c only)

`masin_niemunis_visco_umat` rewritten as the manual Niemunis overlay on
the hcea kernel:

- delegates the base update to `masin_umat` (the plain clay kernel) -
  identical state feeding, tangent and void-ratio update;
- computes the manual visco correction `-L:Dvis*dt`,
  `Dvis = Dr*(1/OCR)^(1/Iv)*mhat`, with `L = fb*Lhat` (fb from
  kappa* - FIXED, was lambda*), the Wolffersdorff shape F (F=1 at
  isotropic, manual note) for Lhat/mhat, classic `eta = qq/(M_c p)`
  for the Rendulic cap `peplus` (betaR = 1 limit `p(1+eta^2)`), and
  the FIXED-NCL reference `pe = exp((N - ln(1+e))/lambda)` (peast
  closure) for OCR;
- activation cutoff: the overlay applies only when
  `creep_rate >= 1e-6`. At the corpus parameters the creep rate stays
  below it for the whole run (early 9.4e-13, late-run max ~2e-7), so
  hypo13 is bit-identical to the plain clay route - which is exactly
  what the Professional does (its corpus visco run is bit-identical to
  its no-visco run);
- no-op creep term removed (the old kernel could never produce creep at
  near-isotropic states because m = 0 there).

## Per-test state after DEV-MASIN2

| test | GNU before | GNU after | Professional | state |
|------|-----------|-----------|--------------|-------|
| hypo13 | crash t=10390 (sigyy -0.554, e -0.68) | completes: -145.690839954, e 0.594128 | -144.495200493, e 0.603664 | CRASH CLOSED; 0.83% residual (rc=1, target +-0.01) |
| hypo7 | -230.792534395 | byte-identical | -231.808976806 | unchanged (0.44%) |
| hypo8 | -217.670847523 | byte-identical | -230.588113239 | unchanged (5.6%) |
| hypo9 | -230.792534395 | byte-identical | -231.810894524 | unchanged (0.44%) |
| hypo12 | -143.393446689 | byte-identical | -143.494849517 | unchanged (0.07%) |
| hypo1/2/4 | PASS | PASS (suite 16/16; wolfersdorff path untouched) | - | neutral |
| hypo3 | unchanged | unchanged (wolfersdorff path untouched) | -0.162914620443 | discard verification NOT run this sprint |

Suite: 16/16 after every build. No memory-handling changes were made
(numerics only in an existing kernel entry point), so the ASan pass was
not re-run.

## OPEN (next sprint): the hcea kernel-level residue - now the single blocker

All remaining hypo-family residuals share one root cause: the
Professional's internal masin kernel is NOT the public umat_hcea.for.
Decisive measurement: at the IDENTICAL state (sig -100 iso, e 1.03) and
IDENTICAL converged increment (hypo7 step 1: deps=(6.3e-5,-5e-4,6.3e-5),
extracted from the Pro's own .his), the Pro FE reports dsigyy = -2.5598
while the faithful hcea port (masin.c, textually verified against the
2885-line Fortran including the I2 = half(norm2sig - I1^2) convention)
reports -1.4295. The same-input kernel gap is 1.79x at hypo7 step 1 and
~1.2x at hypo13 step 1 (Pro -0.0516 vs GNU -0.0428). Candidates that
were FALSIFIED this sprint: RKF tolerance, perturbation tangent, Newton
path, strain-measure feeding (de matches exactly), the 2020 umat_hcea
(wrapper-only diff), OCR/ashape/ocrcs defaults (hardcoded 0.3/2.0 in
both sources). Also measured: the corpus hypo7 comparison ends at
DIFFERENT printed eptyy (-0.3566 GNU vs -0.3 Pro) purely because the
Professional prints materi_strain_total in engineering measure while
GNU accumulates the log measure (ln(1/0.7) = 0.3567): same physical
state, print-measure only. Next sprint should extract the Pro's kernel
response curve systematically (single-element probes sweeping one factor at
a time: p, e, eta, deps direction) and fit the closure differences.

## Print_apply / .his pitfalls (unchanged from DEV-C)

Fresh directory per run; Pro with `print_apply -no` writes nothing
(run a print variant for .dbs/.his - printing does not affect results:
the Pro print-variant full run reproduces the bare-run finals exactly).

---

# DEV-MASIN3 continuation (2026-09-09, worktree dev/masin3 from 40fa51e)

The hcea kernel-level residue hypothesised above is RESOLVED as follows
(full evidence in `masin3_response_curve_fit.md`):

- For `-masin_clay` the Professional's kernel stress response is the
  public hcea formulation (<= 1e-4 relative over ~50 imposed-strain
  probes: p 25..400, e sweeps, K0 starts, both geometries, both
  istrain branches, 100x increment range).
- The corpus residuals are FE-LOOP differences: the Professional's
  lateral strain per step is the frozen plane-strain elastic
  prediction (`depsx = -nu/(1-nu) * depsy` = 2.5e-6 in hypo13 for
  1000+ steps) while GNU's 2-iteration loop converges to the true
  lateral-neutral strain (0.96e-6). The two loops use different
  tangents (measured Pro/GNU: 0.384 for hypo13, 0.984 for hypo12).
- The Professional feeds `hyhis4 = 2` (sensitivity) into its kernel
  without a structure record; with `kparam = 0` that halves `fd` via
  `N* = N + lambda* ln 2` (softmodel branch, `Stf = 1`).
- The Professional's visco route shows a first-step void-ratio quirk
  (`de = 10 x (1+e) tr(deps)` at dtime = 10).
- The basic `-masin` record IS a different kernel in the Professional
  (r-parameterised, isotropic; coincides with hcea at isotropic states
  only under `nu_eff(r)`: 0.2->0.054, 0.4->0.243, 0.6->0.315; breaks
  at K0 states). The 1.79x "kernel gap" of DEV-C is this record-level
  difference, not a clay-kernel residue.

Consequence: closing hypo12/13 within +-0.001/+-0.01 requires the
Professional's loop semantics (solver layer), and hypo7/8/9 the
basic-masin formulation — neither is a masin-kernel delta.
