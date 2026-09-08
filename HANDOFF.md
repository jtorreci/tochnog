# HANDOFF — DEV-B: displacement-primary transient for materi_dynamic<1

Sprint: 4-dev parallel convergence (dynamic1/2/5/8). Branch `dev/dyn` from 19e995a.
Design: ProjectDocs/DIAG-SOLVE-MIXTO.md §17 (measured Professional f=0 map).
Allowed surfaces: calcul.cc, so.cc, general.cc + shared append-only (database.cc,
tochnog.h, tochnog-mod.h) + HANDOFF.md + ProjectDocs/DIAG-SOLVE-MIXTO.md (append only)
+ ProjectDocs/manual-developer/dynamic_displacement_primary.md.

## State (update at every milestone)

- [x] M0 Worktree verified clean @19e995a, branch dev/dyn. HANDOFF.md created.
- [ ] M1 Build green (build_safe 16/16), baseline GNU runs of targets captured.
- [ ] M2 Professional probes: α(s) dt-sweep re-derived; map variants discriminated
      on multi-DOF probes; dynamic5/8 series replication check.
- [ ] M3 Implementation in allowed surfaces, targets green.
- [ ] M4 Blast radius: ground*/undrained/large/matrix families identical.
- [ ] M5 Docs (§17 append + manual-developer/dynamic_displacement_primary.md),
      Engram observation, final commit.

## Key measured facts (to re-verify, source §17 + /tmp/run_dynamic*.log)

- Pro f=0 map (1-DOF, exact to 1e-12):
  u_n = u_{n-1} + dt·v_n; v_n = v_{n-1} + dt·a_n; a_n = α(s)·(F − k·u_{n-1})/m.
- Measured α points: s=0.1→0.999902, 1→0.75, 2.5→0.2568, 4→0.1142, 6→0.0533.
- Candidate closed form (16-pt fit, RE-DERIVE): α(s) = (1+2s²)/(1+s²)².
  Spot check: s=1 → 3/4 ✓; s=4 → 33/289=0.11418 ✓; s=6 → 73/1369=0.053324 ✓.
- Unconditional stability: α·s² → 2 as s→∞ (ω_eff·dt → √2). det(amp)=1, trace=2−αs².
- Targets: dynamic1/2 u(100)=0.095889±1e-4 (17.7 periods → α must match to ~1e-6 rel);
  dynamic5 −8.58577e-5±1e-7; dynamic8 −0.000998525±1e-7.
- GNU failure mode (prior): velocity-staggered, mass on deformed geometry,
  explicit stability limit; dynamic1 drift 6.44767; dynamic5 −8.47308e-5; dynamic8 blowup.

## Evidence paths

- Prior logs: /tmp/run_dynamic1.log ... run_dynamic8.log, /tmp/dyn.log
- Scratch (prior): /tmp/opencode/dyn_pro*, dyn_probes, dyn5
- Pro binary: /tmp/opencode/tn_prof/prof25/linux_64_bit_user_supplied/tochnog_user_supplied
  (LD_LIBRARY_PATH to its dir)
- Targets: /tmp/opencode/tn_prof/tochnog_version_02-08-2026/test/other/dynamic{1,2,5,8}.dat

## Next step

M1: verify build infra + build_safe; baseline GNU target runs.
