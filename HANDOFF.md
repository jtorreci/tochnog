# HANDOFF — DEV-B: displacement-primary transient for materi_dynamic<1

Sprint: 4-dev parallel convergence (dynamic1/2/5/8). Branch `dev/dyn` from 19e995a.
Design: ProjectDocs/DIAG-SOLVE-MIXTO.md §17 (measured Professional f=0 map).
Allowed surfaces: calcul.cc, so.cc, general.cc + shared append-only (database.cc,
tochnog.h, tochnog-mod.h) + HANDOFF.md + ProjectDocs/DIAG-SOLVE-MIXTO.md (append only)
+ ProjectDocs/manual-developer/dynamic_displacement_primary.md.

## State (update at every milestone)

- [x] M0 Worktree verified clean @19e995a, branch dev/dyn. HANDOFF.md created.
- [x] M1 Build green (build_safe 16/16, rc=0, 4.9s incremental).
      Needed: ln -s <mainrepo>/validation-suite validation-suite (gitignored).
      GNU binary needs LD_LIBRARY_PATH=<wt>/external-downloads/numlib-runtime.
- [x] M2 Professional scheme REVERSE-ENGINEERED (machine precision):
      baseline GNU re-verified: dyn1 6.44767, dyn2 0.582799, dyn5 -8.47308e-5,
      dyn8 704490 (all rc=1; target check exits even w/o check_target).
      Professional probes in /tmp/tn_dyn_pro/sweep (1-dof dt sweep + 2-elem bar
      x 3 load patterns). SEE §SCHEME below.
- [ ] M3 Implementation in allowed surfaces, targets green.
- [ ] M4 Blast radius: ground*/undrained/large/matrix families identical.
- [ ] M5 Docs (§17 append + manual-developer/dynamic_displacement_primary.md),
      Engram observation, final commit.

## THE SCHEME (definitive, machine-verified vs Professional binary)

Prior findings PARTIALLY WRONG: doc §17.1 map (a_n = alpha(F-ku_{n-1})/m, u_n =
u_{n-1}+dt v_n) is EXACT only for single-dof; the general scheme is:

    r  = F_n - B^T sigma_t            (sigma frozen at start of step, f=0)
    w1 = (M + dt^2 K)^-1 r
    w2 = (M + dt^2 K)^-1 (M w1)
    a  = 2 w1 - w2
    v += dt a
    u += dt v

with M = lumped mass (row-sum), K = linear elastic stiffness (reference config).
Equivalent: a = alpha(dt^2 M^-1 K) M^-1 r with alpha(x) = (1+2x)/(1+x)^2
(the prior 16-point closed form CONFIRMED EXACT as the 1-dof case: alpha(s^2)
with s^2 = dt^2 lambda of the generalized eigenproblem K phi = lambda M phi).
alpha = 2 beta - beta^2 with beta = (I + dt^2 M^-1 K)^-1 -> the two solves.
Eigenvalues of a_op M match alpha(dt^2 lambda_m) to 2e-5 (fit precision);
commutator [a_op M, M^-1 K] = 3e-14 (machine zero); full-series simulation
matches ALL probes at 5e-13 for dt in {0.05,0.1,0.2,0.5,1.0}, 3 load patterns.
Stability: alpha in (0,1], alpha->1/2 as x->inf, spectral radius = 1
(det=1 exactly): unconditionally stable, zero dissipation.
node_dof print columns (1D, records dis/vel/acc/stress): [disx velx accx sigxx...].
Prior /tmp/run_dynamic2.log was a GNU run (not Pro) - it poisoned the s=0.1 alpha
point (0.990099 is GNU beta=1 map; Pro gives 0.99990197 = closed form).
Evidence: /tmp/tn_dyn_pro/sweep/ (analyze/pooled/hyp/gamma_scan/mapfit/wq/
modaltest .py + probe dirs), /tmp/tn_dyn_base/ (GNU baselines).

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
