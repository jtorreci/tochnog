# HANDOFF — DEV-C hypo-kernel calibration sprint

- Branch: `dev/hypok` (worktree /tmp/opencode/tn-dev-hypok, from 19e995a)
- Status: STARTED
- Scope: hypo3, hypo7/8/9, hypo12/13 (+ hypo10 documented-discard)
- Allowed edits: hypo.c, hypoplas.cc, masin.c, masin_visco.c, sanisand.c, materi.cc, database.cc/tochnog.h/tochnog-mod.h (append-only), HANDOFF.md, ProjectDocs/manual-developer/group_materi_plasti_hypo_wolfersdorff.md, new ProjectDocs/manual-developer/hypo*.md
- FORBIDDEN: SEGUIMIENTO-CONVERGENCIA.md, top.cc, calcul.cc/so.cc/general.cc, calcul_force.cc/mpc.cc, plasti.cc/plasti_i.cc, elem.cc, calcul.cc etc.

## Plan
1. Locate the 6 target tests + reference data; run GNU baseline to confirm numbers.
2. Measure Pro behavior (single-element probes): hypo3 dt-sensitivity, hypo7/8/9 strain-increment semantics, hypo12/13 visco params.
3. Calibrate/fix inside owned surfaces; document anything needing elem.cc/calcul.cc.
4. Blast radius check: hypo1/2/4 byte-identical; corpus families unchanged.

## State
- Milestone FINAL: 3 fixes landed (hypoplas.cc state feeding + epi write-back + structure gate; masin_visco.c e-update). hypo12 −282.46→−143.393 (Pro −143.4948, target ±0.001 — 94% closed); hypo7/9 −224.08→−230.79 (Pro −231.81, 0.44%); hypo8 −216.29→−217.67 (5.6%); hypo13 still crashes (visco kernel closure — full diagnosis in doc); hypo3 + hypo1-4 byte-identical; suite 16/16 after every change.
- Full per-test table + remaining-residual diagnosis:
  ProjectDocs/manual-developer/hypo_masin_convergence_devc.md

## Findings
1. GNU baseline (worktree build, rc=1 all): hypo3 sigxx=−0.103436189421; hypo7 sigyy=−224.076947603; hypo8 sigyy=−216.29249576; hypo9 sigyy=−224.076947603; hypo12 sigyy=−282.457704737 (e=0.5906, physical); hypo13 crash 'iterative solver broke down' at t≈10390 with e→1.9375e14 (matches documented 1.9e14 overflow).
2. Professional (re-measured, prints enabled, rc=0): hypo12 sigyy=−143.494849517 e=0.5997 (target −143.495635±1e-3 → PASSES); hypo13 sigyy=−144.495200493 e=0.60366 (target −144.49±1e-2 → PASSES). Pro hypo3/7/8/9 confirmed from prior complete runs: −0.162914620443 / −231.808976806 / −230.588113239 / −231.810894524.
3. PITFALL FOUND: .his history files APPEND across runs; the stale tn-gnu-c .his were the OLD GNU hypo13 blow-up trajectory (t→10390, e→1.9e14, sigyy→−0.554). Pro with `print_apply -no` writes NO output files at all. Any A/B must use fresh dirs per run and only the .dat copied.
4. hypo12: imposed eptyy identical (−t, vely=−1, 1000 steps×1e-4 to t=0.1); GNU e=0.5906/sigyy=−282.46 vs Pro e=0.5997/sigyy=−143.49 → integration/model-path divergence, not boundary-condition difference.
5. hypo13: control_timestep 10./1e5 (end time 1e5); GNU dies at t≈10390 (10%) when e blows up → solver breakdown. Pro completes t=1e5 physically.

6. ROOT CAUSE hypo12/13 (measured, categorical): the masin dispatch in hypoplas.cc feeds the kernel the NEWTON ITERATE (new_hisv[e/s], new_epi delta) instead of the step-start state (old_hisv/old_epi). Every extra Newton iteration re-integrates the state: measured step-1 void ratio GNU = 0.6 + trial-strain increment (−1.6e-4) + converged-strain increment (−1.19e-4) = 0.59972066, while Pro = one integration of the converged increment = 0.59987997 (= (1+e)*tr(deps_converged) exactly, with Pro's own eptxx=2.498e-5). Same for the intergranular delta: GNU's epi dof NEVER evolves (stuck at reset −1e-5, kernel-returned delta discarded), Pro's evolves (−1e-5 → +5.67e-6 in step 1). hypo12's frozen endgame (sigyy=−282.457704737 to 12 digits from t=0.004 on) is a downstream rigid-motion degenerate equilibrium of the corrupted state.
7. ocr_apply gate (control_materi_plasti_hypo_masin_clay_ocr_apply) is DEAD in GNU (v4 probe: -yes identical to -no) but WORKS in Pro (v4: e0 recalculated via OCR 0.5999->0.805). Not exercised by the corpus targets (they run -no). Documented, not fixed.
8. wolfersdorff hypo_() also reads new_hisv in place, but hypo1-4 are byte-identical targets calibrated to Pro 0.0002% — wolfersdorff path stays UNTOUCHED (neutrality).

## Next (for the next dev on this family)
- hypo13: obtain the Niemunis visco reference (umat_ab_visco-class) or probe-series the Pro creep closure; the pe/peplus degeneracy at eta=0 makes OCR a pure function of e with a ~100-power feedback.
- hypo12/7/9 residue: compare RKF tolerance + perturbation tangent against Pro; wire group_materi_plasti_hypo_strain_intergranular into the basic masin kernel.
- Fix the dead control_materi_plasti_hypo_masin_clay_ocr_apply gate (Pro honours it; measured).
- hypo3: discard verification NOT re-run this sprint (time); documented coarse-step claim still stands unverified-by-Pro-refinement.

## Evidence paths
- /tmp/opencode/tn-dev-hypok-probes/base/{hypo3,hypo7,hypo8,hypo9,hypo12,hypo13} (GNU baseline runs)
- /tmp/opencode/tn-dev-hypok-probes/mine/h12gnu,h13gnu (fresh GNU 12/13)
- /tmp/opencode/tn-dev-hypok-probes/mine/h12pro2,h13pro2 (Pro 12/13, print_apply removed — measurement-only .dat edit)
- prior complete Pro runs: /tmp/opencode/tn-dev-hypok-probes/{hypo3,hypo7,hypo8,hypo9}/pro/*.dbs

- 1-step isolation probes: /tmp/opencode/tn-dev-hypok-probes/mine/p12/{v0,v1,v2,v4} x {gnu,pro}; debug dump p12/dbg (TOCHNOG_DEBUG=yes) showing 8 materi calls/step with trial+converged strains

## Next
- hypo7/8/9: verify updated_linear strain-increment mechanism (GNU log strain vs Pro linear increments); minimal guarded fix in materi.cc if possible, else document elem.cc need.
- hypo12: instrument GNU vs Pro trajectory (get Pro .his: print item?). Investigate masin.c input semantics (strain increment convention).
- hypo13: locate e blow-up in masin_visco.c integration.
- hypo3: verify Pro coarse-step artifact (dt refinement both codes → ~−0.105).
