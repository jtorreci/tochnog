# Hypoplasticity kernel calibration (2026-09-07, DEV-C pilot)

Reference measured with the Professional binary 25-10-2023
(`tn_prof/prof25/.../tochnog_user_supplied`) on the corpus
`test/other/hypo*.dat`; every number below is read from the final
`node_dof`/`post_point_dof` records of the `.dbs` written next to the
`.dat` (full precision) and, for the step series, from truncation runs
(the same `.dat` with the `control_timestep` end time shortened, one run
per sample).

## Code changes applied (commit `feat(hypo)`)

`hypo.c` only (no new keywords; `database.cc`/`tochnog.h` untouched):

1. **Intergranular strain: clamp `rho = ||S||/R` to `<= 1`.**
   - Justification: Professional manual, intergranular strains extension:
     "From the evolution equation it follows that rho must remain between
     0 and 1." The discrete midpoint rule overshoots `rho > 1` whenever the
     strain direction rotates inside a step (e.g. the isotropic reset
     `S = (-1,-1,-1) R/sqrt(3)` of hypo2 in an oedometric path). Unclamped,
     `rho > 1` gives a NEGATIVE weight to `mR` in the stiffness
     `[rho^chi*mT + (1-rho^chi)*mR]`.
   - Effect: hypo2 `sigxx` -0.176 -> -0.17447 (Professional -0.17284),
     still rc=0. hypo1/3/4 unaffected (hypo3's reset is aligned with the
     loading direction, so rho never overshoots there).

2. **Substep tolerance of the stress-increment criterion: 1% -> 2%.**
   - The kernel halves its substep `dt` while the trial substep stress
     increment exceeds a fraction of the current stress magnitude. The
     historical constant (hypo.f, GNU 2014) was 1%. Measured against the
     Professional binary the *effective* tolerance of the Professional
     kernel is ~2%: hypo1 `sigyy` -862.766 @1% vs Professional -862.929;
     @2% the GNU gives -862.927 (target -862.92 +- 0.1 -> rc=0); @5% it
     overshoots to -863.29. The sweep in [1%, 2.8%] keeps hypo1 rc=0.

## Results per corpus test (GNU before -> after, vs Professional)

| test | dof | GNU before | GNU after | Professional | rc | change |
|------|-----|-----------|-----------|-------------|----|--------|
| hypo1 (wolfersdorff) | sigyy | -862.766345621 | -862.927047891 | -862.9290976 | 1 -> **0** | substep tol 2% |
| hypo2 (wolfersdorff+epi) | sigxx | -0.176 | -0.174467815 | -0.172844871 | 0 | rho clamp (in tol both) |
| hypo3 (wolfersdorff+epi, K0 unload) | sigxx | -0.104030562 | -0.103436189 | -0.162914620 | 1 | unchanged (see below) |
| hypo4 (wolfersdorff+epi cyclic) | sigxy | 0.010526271 | 0.010526271 | 0.010505331 | 0 | unchanged |
| hypo7 (masin) | sigyy | -224.0769 | -224.0769 | -231.8089 | 1 | unchanged (see below) |
| hypo8 (masin+structure) | sigyy | -216.2925 | -216.2925 | -230.5881 | 1 | unchanged |
| hypo9 (masin+ig) | sigyy | -224.0769 | -224.0769 | -231.8109 | 1 | unchanged |
| hypo12 (masin clay visco Dr/Iv) | sigyy | -282.458 | -282.458 | -143.4956 | 1 | unchanged (see below) |
| hypo13 (Niemunis visco) | sigyy | crash (e -> 1.9e14) | same | -144.49 | 1 | unchanged (see below) |

Blast radius re-verified rc=0: hypo1, hypo2, hypo4, mohr_coul_direct1-3,
validation_12, visc_pl1. ground11_phreatic_level / ground11_nonsaturated
(Lapack direct solver singular) and direct_shear_drained_hypo (unknown
keyword `incremental_driver`) fail identically on the pristine kernel
(pre-existing, not hypo.c related).

## Diagnosed as scheme/upstream gaps (no kernel change)

### hypo3 — reversal step integrated at 10*R per step (36% deviation)

Truncation series show the loading branch (t = 0 -> 0.01) matches the
Professional within 0.4% (sigxx -0.33446 vs -0.33592 at t=0.01); the whole
36% gap opens in the TWO unloading steps (t=0.010 -> 0.012). The reversal
step carries a strain increment of 1e-3 = **10 times the intergranular
radius R** (R=1e-4): the epi flips from -R to +R inside a single step, a
marginally resolved integration.

Evidence that this is NOT a kernel-equation gap:
- Running both codes with finer FE steps (dt=1e-4, 5e-5) converges them to
  the SAME value: GNU -0.1037/-0.1040, Professional -0.1084/-0.1063
  (converged solution of the model ~= -0.105). The Professional's corpus
  value -0.16291 is its own COARSE-step result (it integrates the reversal
  with 1-2 internal substeps vs ~16 for the GNU 1% criterion).
- Without intergranular strain (same .dat, epi machinery removed) both
  codes agree to 0.8% (GNU -0.1685, Professional -0.1672): the base
  wolfersdorff branch is consistent.
- Sweeping the GNU substep tolerance from 1% to 100% never reproduces
  -0.163 (values plateau at -0.098/-0.154 and diverge beyond ~40%): the
  Professional's coarse-step reversal response is not reachable by this
  tolerance knob without making the reversal integration unstable.

Conclusion: the corpus target (-0.1618, Herle special-purpose program) and
the Professional value (-0.16291) are coarse-step artifacts of the reversal
integration; the converged model value of both kernels is ~-0.105. hypo4
(cyclic shear, strain per step = 0.1*R) agrees between the codes to 0.2%,
confirming the models match whenever the steps resolve R. A fix would
require reproducing the Professional's coarse-step epi-flip integration
(scheme-level, not a parameter).

### hypo7/8/9 — masin basic/clay tests under -updated_linear (3-6%)

The strain state accumulated by the FE differs between the codes in the
SAME .dat (axisymmetric + materi_velocity_integrated +
group_materi_memory -updated_linear, hypo7):
- Professional: eptyy = -0.300000000000 (exactly linear in time; each step
  -5e-4), eptxx = +0.16950.
- GNU: eptyy = -0.35657 (= -ln(1-0.3): the strain increments are evaluated
  in the CURRENT geometry each step, accumulating logarithmic strain),
  eptxx = +0.19004.

The void ratio evolved by the masin kernel is consistent with the LINEAR
strain path in the Professional (e: 1.03 -> 1.1107, de ~= (1+e)*tr(ept))
and with a different (larger) volumetric path in the GNU (e: 1.03 ->
1.1263, dip to 0.941 vs Professional dip to 1.006). The sigyy divergence
(3-6%) is a downstream consequence of the different imposed strain path,
not of the masin kernel: removing materi_strain_plasti changes nothing;
the masin.c kernel itself is a validated port of umat_hcea.for (5e-7 at
driver level). The -updated_linear kinematics/velocity-integrated strain
accumulation lives upstream of the hypoplasticity dispatch (not in
hypo.c/hypoplas.cc/masin.c) and must be aligned there.

### hypo12/hypo13 — Masin clay visco (Niemunis Dr/Iv law)

hypo12 runs to completion (sigyy -282.458, e 0.5906 stable) but is 2x off
the Professional target -143.4956; hypo13 diverges the void ratio (e ->
1.9e14, "severe error in Masin hypoplasticity"). Both use the Niemunis
visco law (group_materi_plasti_hypo_masin_clay_visco, Dr=1e-6 Iv=0.1),
which is implemented from the manual theory only (no reference UMAT
exists; see group_materi_plasti_hypo_masin.md). Calibrating it would
require a numeric reference for the Dr/Iv creep formulation; per the pilot
rule the divergence is documented, not forced.
