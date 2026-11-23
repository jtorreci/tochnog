# materi_dynamic / control_materi_dynamic

Developer notes for the dynamics blending factor.

## Where

- `tochnog.h` / `tochnog-mod.h`: `MATERI_DYNAMIC` record enum (global, no
  index) next to `CONTROL_MATERI_DYNAMIC`.
- `database.cc` `db_initialize()`: `materi_dynamic` registered as a
  DOUBLE_PRECISION global record (`no_index = 1`, class CONTROL).
  `control_materi_dynamic` was already registered but as INTEGER; changed to
  DOUBLE_PRECISION because the factor is a real number (e.g. 0.0 in
  earthquake_2).
- `materi.cc` (continuum elements): after the constitutive stress and the
  SRI adjustments, `sigvec` (the momentum internal-force stress) is blended
  with the integration-point stress of the previous converged step
  (`old_sig`): `sigvec = (1-f)*old_sig + f*sigvec`. The momentum stiffness
  (`stiffness`, SRI `stiffness_shear`, and the SRI shear feedback force) is
  scaled by `f`. The stress *state* dofs (`new_sig`, line ~860) keep the
  full update, and `options_element_dof` writes stay untouched.
- `truss.cc` (truss/bar elements): after the force caps (plasti/rope), the
  momentum force is `new_truss_force_mix = old_truss_force +
  f*(new_truss_force - old_truss_force)` and `truss_stiffness` is scaled by
  `f` for the matrix. `element_truss_force` (VERSION_NEW state) keeps the
  full `new_truss_force`.
- Factor resolution in both routines: global `materi_dynamic` first, then
  `control_materi_dynamic icontrol` overrides when the current control
  index has the record. Out-of-range factors (< 0 or > 1) -> `db_error`.

## Semantics

Per the manual (6.800): `sigma_used = (1-f)*sigma(t) + f*sigma(t+dt)` with
default f = 1 (fully implicit). Implemented incrementally in the momentum
equation, this is the generalized trapezoidal family parameterized by f:
f = 1 reproduces the historic GNU scheme byte-for-byte (verified: the
recovered-suite checks and truss5 stay rc=0, and the f = 1 code path only
adds record reads).

## Measured status (Professional binary 25-10-2023 A/B)

| test | GNU f=1 (implicit) | GNU f=0 | Professional f=0 | Professional f=1 |
|------|--------------------|---------|------------------|------------------|
| dynamic1 (truss) | 1.912 | 6.45 (growing) | 0.095889 (target) | 0.9955 |
| dynamic2 (bar2) | - | 0.582799 | 0.095889 (target) | - |
| dynamic5 (quad4) | -7.29e-5 | -8.47e-5 | -8.58577e-5 (target) | - |
| dynamic8 (quad4) | -3.48e-4 | diverges (704490) | -9.985e-4 (target) | - |

Conclusions (2026-09-07 fine diagnosis, DIAG-SOLVE-MIXTO S17 — NO code
change; the numbers below are the full measured mechanism):
1. f = 0 moves the GNU response toward the Professional targets (dynamic5:
   15% -> 1.3% error), confirming the record semantics (frozen sigma_t in
   the momentum, stiffness*f in the matrix).
2. The Professional f=0 is a displacement-based semi-implicit map
   (measured to 1e-13 on the state series of dynamic1/2): u_n =
   u_{n-1}+dt*v_n, v_n = v_{n-1}+dt*a_n, a_n = alpha(s)*(F - k*u_{n-1})/m
   with a step-dependent relaxation alpha(s) ~ 1 (s=w*dt << 1), alpha(s)
   -> 2/s^2 (large s), i.e. w_eff*dt -> sqrt(2): UNCONDITIONALLY stable
   (the Pro runs dynamic2 at dt=2.5/4/6 bounded and dynamic8 at 20x the
   corpus dt rc=0). Mass is constant (reference configuration); the load
   "bounda_time <icontrol> 1." with ONE value = CONSTANT magnitude 1 for
   the whole run (not a pulse).
3. The GNU f=0 is the velocity-staggered explicit map (matrix = M/dt
   only), with TWO equilibrium iterations per step and the continuum mass
   integrated on the iterate-deformed geometry (assembled diag 10 ->
   10.1 -> 10.198 as the bar stretches). Differences vs the Pro:
   (a) explicit stability limit w*dt < 2 (GNU dynamic2 at dt=2.5
   diverges to -2678 where the Pro is bounded; dynamic8 diverges);
   (b) updated kinematics with no finite equilibrium at the corpus loads
   F = EA (dynamic1/2 drift); (c) small residual dissipation/phase even
   at small strain (dynamic5: amplitude 0.8% low; dt-refinement is
   first-order and the dt->0 extrapolation -8.52e-5 still misses the
   target -8.58577e-5 +-1e-7).
4. The corpus targets (±1e-4 at t=100 for dynamic1/2, ±1e-7 for
   dynamic5/8) pin the Professional's DISCRETE trajectory: dynamic1/2
   would need the integrator frequency matched to ~1e-6 relative (even a
   nominal symplectic replica misses by 0.0021 at t=100 because the Pro's
   alpha(0.1) = 0.999902 != 1). A consistent closure inside the GNU solve
   (the mpc3/4 precedent) does NOT transfer: the gap is the transient
   FORMULATION (reference-configuration mass/stretch + semi-implicit
   displacement solve), i.e. a displacement-primary time integration
   branch — scheme-level work unit with a wide blast radius.
5. The dynamics family (dynamic1/2/5/8) stays RUNFAIL with this measured
   root cause: SCHEME/FORMULATION (see DIAG-SOLVE-MIXTO S17), NOT parse
   and NOT the record semantics.

## Pending

- Apparent-modulus based damping `support_edge_normal_damping_automatic`
  variants read the material state; unaffected by the blend.
- Full displacement-primary transient branch for materi_dynamic < 1 with
  the Pro's measured map (scope in DIAG-SOLVE-MIXTO S17.4).
- SMC reading: implemented in `bounda.cc` (`bounda_time_smc_read`), see bounda_time_smc.md.

## Effective factor helper (F4 T6d / C11, 2026-10-11)

- `general.cc` `materi_dynamic_effective( long int *specified )`: the
  single source of the factor. Global `MATERI_DYNAMIC` (default 1), then
  `CONTROL_MATERI_DYNAMIC` at `ICONTROL` overrules it; `*specified` = 1
  when one of the two records applies. Returns 1 while
  `disp_primary_reference == -YES` (the f=1 auxiliary reassembly of the
  displacement-primary scheme). Declared in `tochnog.h` / `tochnog-mod.h`.
- Readers: `so.cc` `disp_primary_gate()` (requires `specified` and
  0 <= f < 1; the old exclusion of `CONTROL_MATERI_DYNAMIC` is gone) and
  `solve()` (`disp_primary_f`); `truss.cc` (reference-configuration test
  and the force blend), `materi.cc` (reference test and momentum blend),
  `elem.cc` (reference mass). `solve()` no longer PUTs `MATERI_DYNAMIC = 1`
  around the auxiliary `element_loop()` and restores it: the marker
  `disp_primary_reference` alone selects f = 1, so a control record can no
  longer override the temporary value.
- Verified against the Professional 02-08-2026 (node 2 disx): dynamic1 +
  only `control_materi_dynamic 10 0.` 6.4477 -> 0.0958890; global +
  control 0. 0.0937873 -> 0.0958890; control 0.5 1.45599 -> 0.932818; no
  dynamic record 1.91225 -> 0.995505 (the last one through the truss
  memory default, see `group_truss_memory.md`); dynamic3 32.6583 ->
  60.44025.

## Differences with the Professional version

- None known for the factor resolution: the former C11 difference
  (`control_materi_dynamic` ignored by the displacement-primary gate,
  truss default memory) was resolved on 2026-10-11 (F4 T6d; BACKLOG C11).
  Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
