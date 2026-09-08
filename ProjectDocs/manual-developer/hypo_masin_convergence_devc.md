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
