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
- Milestone 1 (measurement): worktree build OK (suite 16/16 green). GNU baseline reproduced documented numbers. Pro re-measured for all 6 targets.

## Findings
1. GNU baseline (worktree build, rc=1 all): hypo3 sigxx=−0.103436189421; hypo7 sigyy=−224.076947603; hypo8 sigyy=−216.29249576; hypo9 sigyy=−224.076947603; hypo12 sigyy=−282.457704737 (e=0.5906, physical); hypo13 crash 'iterative solver broke down' at t≈10390 with e→1.9375e14 (matches documented 1.9e14 overflow).
2. Professional (re-measured, prints enabled, rc=0): hypo12 sigyy=−143.494849517 e=0.5997 (target −143.495635±1e-3 → PASSES); hypo13 sigyy=−144.495200493 e=0.60366 (target −144.49±1e-2 → PASSES). Pro hypo3/7/8/9 confirmed from prior complete runs: −0.162914620443 / −231.808976806 / −230.588113239 / −231.810894524.
3. PITFALL FOUND: .his history files APPEND across runs; the stale tn-gnu-c .his were the OLD GNU hypo13 blow-up trajectory (t→10390, e→1.9e14, sigyy→−0.554). Pro with `print_apply -no` writes NO output files at all. Any A/B must use fresh dirs per run and only the .dat copied.
4. hypo12: imposed eptyy identical (−t, vely=−1, 1000 steps×1e-4 to t=0.1); GNU e=0.5906/sigyy=−282.46 vs Pro e=0.5997/sigyy=−143.49 → integration/model-path divergence, not boundary-condition difference.
5. hypo13: control_timestep 10./1e5 (end time 1e5); GNU dies at t≈10390 (10%) when e blows up → solver breakdown. Pro completes t=1e5 physically.

## Evidence paths
- /tmp/opencode/tn-dev-hypok-probes/base/{hypo3,hypo7,hypo8,hypo9,hypo12,hypo13} (GNU baseline runs)
- /tmp/opencode/tn-dev-hypok-probes/mine/h12gnu,h13gnu (fresh GNU 12/13)
- /tmp/opencode/tn-dev-hypok-probes/mine/h12pro2,h13pro2 (Pro 12/13, print_apply removed — measurement-only .dat edit)
- prior complete Pro runs: /tmp/opencode/tn-dev-hypok-probes/{hypo3,hypo7,hypo8,hypo9}/pro/*.dbs

## Next
- hypo7/8/9: verify updated_linear strain-increment mechanism (GNU log strain vs Pro linear increments); minimal guarded fix in materi.cc if possible, else document elem.cc need.
- hypo12: instrument GNU vs Pro trajectory (get Pro .his: print item?). Investigate masin.c input semantics (strain increment convention).
- hypo13: locate e blow-up in masin_visco.c integration.
- hypo3: verify Pro coarse-step artifact (dt refinement both codes → ~−0.105).
