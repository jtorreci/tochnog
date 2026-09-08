# DEV-A HANDOFF — incremental_driver

Sprint: 4-dev parallel convergence (GNU tochnog vs Tochnog Professional 25-10-2023).
Branch: `dev/incdrv` (from 19e995a). Worktree: /tmp/opencode/tn-dev-incdrv.

## State
- [x] Worktree verified clean @ 19e995a.
- [x] HANDOFF created first.
- [ ] Read DECISIONES.md + SEGUIMIENTO notes on incremental_driver.
- [ ] Re-verify prior findings from /tmp/incdrv_* artifacts.
- [ ] Finish Pro algorithm enumeration → precise spec in incremental_driver doc.
- [ ] Pro A/B probes before implementation.
- [ ] Implement in incremental.cc/top.cc/database.cc/tochnog.h/tochnog-mod.h.
- [ ] Suite 16/16 green; corpus targets PASS.

## Targets (test/incremental_driver/ of Pro test tree)
- incremental_driver_syntax.dat — currently PARSE
- oedometric_drained_mohr_coulomb.dat — RUNFAIL
- triaxial_compression_drained_mohr_coulomb.dat — RUNFAIL
- triaxial_compression_undrained_mohr_coulomb.dat — RUNFAIL (needs groundflow coupling)
- NOT ours: direct_shear_drained_hypo.dat, triaxial_compression_drained_hypo_intergranular_isa.dat (DEV-C kernel)

## Key prior findings (RE-VERIFY, from killed attempts; raw: /tmp/incdrv_*)
- Pro generates the driver model INTERNALLY from the .dat declaration (probeta mesh etc.); writes incremental_driver_result.txt.
- Core Pro functions identified from asm: mohr_coul_direct, deften_decompose; undrained needs groundflow coupling.
- Probe result: /tmp/opencode/probe_incdrv/incremental_driver_result.txt

## Evidence paths
- (to be updated per milestone)

## Commits
- (to be appended at every milestone)
