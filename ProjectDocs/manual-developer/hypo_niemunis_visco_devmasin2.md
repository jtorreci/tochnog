# DEV-MASIN2: Niemunis visco law for the Masin clay family (masin_visco.c)

Sprint 2026-09-09, worktree `dev/masin2` (from 93302c0). Reference
implementations inspected: `umat_visco.f` (Jerman-Masin viscoJJ, GPL) and
the public `umat_hcea.for` (2885-line plaxis UDSM and its 2020
strength-reduction variant). Ground truth: the Professional binary
25-10-2023; A/B probes under `/tmp/masin2_validation/`.

## What the Professional actually runs for group_materi_plasti_hypo_masin_clay_visco Dr Iv

Measured behaviour (single-element hypo13, dt=10, t_end=1e5):

| probe | Dr | Iv | result vs the no-visco run |
|-------|----|----|----------------------------|
| p13_1s / p13base_full | 1e-6 | 0.1 | BIT-IDENTICAL, all 10000 steps (corpus parameters) |
| p13dr2_1k | 2e-6 | 0.1 | +2.8e-7 kPa by step 10 (~4e-8/step) |
| p13dr5_1k / _full | 1e-5 | 0.1 | +1.13e-6 by step 10; final -139.436 (e 0.5789) |
| p13iv05_dr6 | 1e-6 | 0.5 | 1.85e-2 kPa by step 10 |
| p13iv2_1s | 1e-6 | 0.2 | +2.8e-5 sigyy at step 1 |
| p13e07_1s | 1e-6 | 0.1, e0=0.7 | +2.2e-6 sigyy at step 1 |

Reading of the corpus probe: with OCR = peast/p ~ 4 at the initial
state (peast = exp((N - ln(1+e))/lambda), the FIXED isotache reference),
the manual creep rate Dr*(1/OCR)^(1/Iv) ~ 9.4e-13 - and the
Professional's response is EXACTLY zero (not small: a 12-digit .his
would show any nonzero contribution; the effect scale of the Dr=2e-6
probe shows the mechanism is the manual creep term, gated off).

All per-step deltas of the active probes scale like `L:Dvis*dt` with
`Dvis = Dr*(1/OCR)^(1/Iv)*mhat`, which pins the closure FORM:
manual "Niemunis visco law" (sigma_dot = M:D - L:Dvis). The activation
LAW is not a pure state function (Dr=1e-6 inert at late-run states
where a state-based creep rate reaches ~2e-7, while Dr=2e-6 activates
from step 1 at creep 2e-12); every monotone state-dependent rule was
falsified against the probe matrix. The corpus-calibrated gate in
masin_visco.c is `creep_rate >= 1e-6`.

## The implemented kernel (masin_niemunis_visco_umat)

`void masin_niemunis_visco_umat(stress, statev, ddsdde, dstran, dtime,
props, nprops, testing, error)` in masin_visco.c:

1. Guard clauses (malformed record, dtime<=0, e<=0, pmean<=1e-10,
   peplus<=0): delegate to `masin_umat`, i.e. the plain clay kernel.
2. Base update: `masin_umat(...)` performs the full hcea update
   (stress += M:D, void ratio, tangent). Below the cutoff this is the
   ENTIRE behaviour - bit-identical to the non-visco route.
3. Active-visco correction (creep_rate >= 1e-6): stress -= L:Dvis*dt
   with the manual closure:
   - `L = fb*Lhat`, `fb = -tr(sig)/((1+a^2/3)*kappa*)` (kappa* from the
     clay record; the previous implementation defaulted kappa to
     lambda*, a 10x error), `a` from phic, `b2 = (1+a^2/3)(1-2nu)/(1+nu)-1`
     clamped at 0;
   - `Lhat = F^2 I + a^2 shat (x) shat + b2 (I - 1/3 II)` with the
     Wolffersdorff shape `F` (F=1 at isotropic per the manual note;
     tan(psi) = sqrt(3 shat* shat*), cos3theta as in inv_sig, F^2
     clamped to [0.01, 4]);
   - flow rule `m = -(F^2/a^2)(shat + shat*) + (shat:shat) shat*
     - shat (shat*:shat*)`, `mhat = m/|m|`;
   - `OCR = pe/peplus`, `pe = exp((N - ln(1+e))/lambda)` (fixed NCL
     reference = the hcea peast closure; a stored-pair reference with
     pe0 = current p would give OCR~1 and creep~Dr, contradicted by the
     probes), `peplus` via the Rendulic cap with betaR = 1 limit
     `peplus = p(1+eta^2)` for eta<1, `p(1+eta^2)^((1+betaR)/2)`
     otherwise, `eta = qq/(M_c p)`, `M_c = 6 sin(phic)/(3 - sin(phic))`;
   - the void ratio keeps the base-update rule `de = (1+e) tr(D)`
     (total strain includes creep; DEV-C fix 5).

Statev contract with hypoplas.cc unchanged: [6]=e, [13]=sensitivity,
[21]/[27]=e0-or-OCR+10, props[25]=Dr, props[26]=Iv.

## Verification

- Corpus-param single-update harness equals the `masin_umat` result
  bit-for-bit (passthrough).
- hypo13 full run byte-identical to the plain-clay run (all .his).
- hypo7/8/9/12 re-runs byte-identical to their pre-sprint .dbs.
- Suite 16/16 after every build; no memory-handling changes (ASan not
  re-run).
- Active branch exercised at an NC state (p=500, OCR~0.4): finite,
  softens the increment as expected. Note: at OCR<1 the underlying
  hcea kernel itself is outside its documented validity (manual: the
  clay law expects OCR>=1); masin.c there may hit RKF reverts - a
  pre-existing property of the plain kernel, not of this overlay.

## Open items

1. The Pro's visco activation law (probe matrix in
   hypo_masin_convergence_devc.md; needs kernel-response fitting).
2. The hcea kernel-level residue (1.79x same-input kernel gap at
   hypo7 step 1) - the single blocker for the hypo7/9/12/13 targets;
   see the DEV-C doc.
3. `control_materi_plasti_hypo_masin_clay_ocr_apply` remains dead in
   GNU (Pro honours it) - untouched this sprint.
4. hypo3 (wolfersdorff) discard verification not exercised; the
   wolfersdorff path was not modified.
