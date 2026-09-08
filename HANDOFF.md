# HANDOFF — DEV-D: force16 / force17 / mpc9 residuals

Sprint: 4-dev parallel convergence (GNU → Professional, binary 25-10-2023).
Branch: dev/fmpc (from 19e995a), worktree /tmp/opencode/tn-dev-fmpc.
Started: 2026-02-06 (sprint slot 45 min corpus window; commands <10 min).

## State

- [x] Worktree verified clean on dev/fmpc @ 19e995a.
- [ ] Read DIAG-SOLVE-MIXTO §16-17 + SEGUIMIENTO dev/segv + dev/veln entries.
- [ ] Locate corpus inputs force16.dat / force17.dat / mpc9.dat.
- [ ] Measure Pro reference values for the three tests.
- [ ] Baseline GNU results (expect RUNFAIL with clear diagnoses).
- [ ] mpc9: re-verify prior measurement (all four nodes at 0.4545 = 1/2.2, condensed diag K22+K33 = 1.833333, v2 = 5/11), implement in mpc.cc only.
- [ ] force16: decode Pro force16 semantics, implement in calcul_force.cc.
- [ ] force17: implement consumption of `_force`/`_inertia` variants (calcul_force.cc; area.cc UNOWNED — document if needed).
- [ ] Blast radius: mpc*/force*/post* families before/after; mpc1/2/7 byte-identical.
- [ ] build_safe suite 16/16 green.
- [ ] Final commit + Engram observation.

## Findings

(to be filled with measurements; evidence paths in /tmp scratch dirs tn-dev-fmpc-*)

## Next step

Read the two ProjectDocs sections, then locate corpus inputs and measure Pro.

## Evidence paths

- Pro binary: /tmp/opencode/tn_prof/prof25/linux_64_bit_user_supplied/tochnog_user_supplied
- Corpus: /tmp/opencode/tn_prof/tochnog_version_02-08-2026/test/
- Scratch: /tmp/tn-dev-fmpc-*
