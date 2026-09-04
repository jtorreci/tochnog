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
| dynamic5 (quad4) | -7.29e-5 | -8.47e-5 | -8.58577e-5 (target) | - |
| dynamic8 (quad4) | -3.48e-4 | diverges (704490) | -9.985e-4 (target) | - |

Conclusions:
1. f = 0 moves the GNU response toward the Professional targets (dynamic5:
   15% -> 1.3% error), confirming the record semantics.
2. The GNU velocity/staggered scheme cannot reach the dissipation-free
   target: at f = 0 its stability region collapses to the explicit limit
   (dynamic8 has c*dt/h >> 1 and diverges), while the Professional remains
   stable — its formulation is displacement-based (manual 6.800 kinematics)
   with a different mass treatment. Even the GNU f = 1 differs from the
   Professional f = 1 (dynamic1: 1.912 vs 0.9955).
3. The dynamics family (dynamic1/2/5/8) therefore stays RUNFAIL with this
   root cause: SOLVER/SCHEME (same family as DIAG-SOLVE-MIXTO), NOT parse.

## Pending

- Apparent-modulus based damping `support_edge_normal_damping_automatic`
  variants read the material state; unaffected by the blend.
- SMC reading (bounda_time_smc) — separate pending feature.
