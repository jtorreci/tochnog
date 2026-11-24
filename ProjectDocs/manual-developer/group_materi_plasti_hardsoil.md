# group_materi_plasti_hardsoil

## Implementación

- **Ley**: new block in `plasti_rule()` in `plasti.cc`, right after the
  `GROUP_MATERI_PLASTI_CAP1` block. It reads the group record
  `GROUP_MATERI_PLASTI_HARDSOIL` (DOUBLE, fixed length 4) and
  implements both `GET_YIELD_RULE` and `GET_FLOW_RULE` (associative:
  `f_flow = f_yield`) for `plasti_type == -NONE` and for the explicit
  type. The flow gradient is provided by the standard driver via central
  finite differences.
- **Keywords** registered in `database.cc`:
  - `group_materi_plasti_hardsoil` (DOUBLE, length 4, required
    GROUP_TYPE): `phi c psi Rf`.
- **check.cc**: `GROUP_MATERI_PLASTI_HARDSOIL` requires
  `materi_stress`, `materi_velocity`, `materi_plasti_kappa` (the
  hardening variable gamma_p), `materi_strain_plasti_hardsoil` and
  `materi_plasti_hardsoil_history` (the manual prescribes both
  initializations).
- The block requires `GROUP_MATERI_ELASTI_HARDSOIL` (error + exit if
  missing: the yield function needs E50/Eur).

## Física

```
f = q/(E50*(1 - q/qa)) - 2*q/Eur - gamma_p
q = sqrt(0.5*[(s11-s22)^2 + (s22-s33)^2 + (s33-s11)^2] + 3*(s12^2+s23^2+s31^2))
qa = qf/Rf
qf = 2*sin(phi)*(sig3 + c*cot(phi))/(1 - sin(phi))     (derived from Mohr-Coulomb at failure)
```

- **sig3**: same mapping as the elastic block (largest algebraic
  eigenvalue = the manual's least-compressive sig3; the base is
  `-(largest eigenvalue) + c*cot(phi)`).
- **E50/Eur**: power laws at the CURRENT trial stress `sig[]` (the
  coupled elasticity-plasticity; one step lag documented on the elastic
  page). Same clamps: base `<= 0` -> E = Eref; `m == 0` -> E = Eref;
  non-positive reference base -> db_error.
- **gamma_p (documented decision)**: the hardening variable is the
  `materi_plasti_kappa` dof (`kappa = int sqrt(0.5*deps_p:deps_p)`, the
  measure the manual calls "equivalent plastic shear strain"; the yield
  function is calibrated to it) PLUS the extra initial contribution of
  `control_materi_plasti_hardsoil_gammap_initial` read from the element
  record `ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL`
  (`GET_IF_EXISTS`, buffer initialized to 0; values `>= 0` after the
  first-timestep initialization, negative sentinel values are treated as
  0 by clamping at the buffer initialization). The manual defines NO
  node dof for gamma_p (its record 4.22 is abs(p)); the kappa coupling is
  the port decision, documented.
- **psi**: accepted but unused (the flow is associative: `f_flow =
  f_yield`); the non-associative HS flow potential is PENDING.
- **Singularity q >= qa**: the surface is asymptotic to `qa` (the stress
  can never reach it); `1 - q/qa <= 0` would give a NEGATIVE f (a trial
  beyond the asymptote would look elastic — wrong). The block sets
  `f = DBL_MAX/1.e6` (the driver's "too large time steps" exit
  threshold, NOT exceeded) so the return is driven hard; if the plastic
  iterations cannot bring q below qa the driver exits with "too large
  time steps". `qa <= 0` (Rf <= 0 or a non-positive base) -> db_error.
  Base `<= 0` -> `f = NO_YIELD_F` (no yielding from the hardsoil law;
  use a tension cut-off for those states, documented).
- **Flow direction**: with `m = 0` the yield depends on `q` only, so the
  numerical gradient is purely deviatoric and the return preserves the
  mean pressure (verified in the tests: the relaxation keeps the mean at
  the initial value).

## Gotchas (numeric, IMPORTANT)

- The HS yield function is O(q/E) ~ 1e-4 (vs O(10-100) for cap1/DP): the
  cutting plane converges at the INTEGRATION POINT level (f -> ~1e-7),
  but the STRESS DOF of the coupled system does NOT robustly propagate
  the returned value (the accumulated dof update converges to the
  elastic trial state). Consequence: a test with continuous loading from
  the virgin state shows the ELASTIC response in the dofs, while the
  relaxation tests (initial deviatoric stress, no loading) DO show the
  return in the dofs (partial relaxation + kappa growth) with the end
  state ON the surface (f ~ 0). The tests validate the surface and the
  hardening through the relaxation A/B family; the loading-path
  validation is PENDING (would need the incremental driver or a
  reformulated stress update).
- `mhardsoil_plast` (initial sigxx = -4): kappa 0.00104991, sigxx
  -3.11898 (after the F4 T5 kappa fix below; the old 0.0006 / -2.44 was
  the even-iteration branch of an oscillation).

## Kappa read inside the return mapping (F4 T5 follow-up, 2026-10-10)

- **Defect**: `plasti.cc` `plasti_rule()` reads `gamma_p` from
  `new_unknowns[kap_indx]` (hardsoil, line ~389). `new_unknowns` is the
  interpolated NODAL kappa that the previous material evaluation of the
  SAME step wrote through `parallel_new_dof_diagonal`. Each evaluation
  therefore saw the plastic strain of the previous one as history: with
  zero total strain `mhardsoil_plast` alternated kap 0.0030 (1 or 3
  equilibrium iterations) / 0.0006 (2, 4, 8, 20, 50). The post-solve
  material state pass (F4 T5) adds one evaluation and flipped the branch.
- **Fix** (`stress.cc` `set_stress()`, plastic iterations): the yield
  rules get `pl_unknowns`, a copy of `new_unknowns` whose `kap_indx` and
  `kapsh_indx` slots hold the step-start value (`old_kappa`,
  `old_kapsh`) and, after every local iteration, `old + increment` of
  this return (the `new_kappa` / `new_kapsh` the driver already
  computes). The four `plasti_rule()` calls of the loop use it (the
  central-difference flow gradient passes it on).
- **Scope**: every yield rule that reads kappa or kappa_shear through
  that array: hardsoil `gamma_p`, the kappa-hardened Mohr-Coulomb
  variants (`plasti.cc` ~774/776 and ~1118) and the von Mises Nadai
  hardening (~1316).
- **Same pattern, NOT changed** (listed in BACKLOG C51): the cap
  `epp_vol` (`new_unknowns[epp_indx]`, ~221), the cap1 `pc`
  (`new_unknowns[cap1_indx]`, ~282), the void fraction (~957/987) and
  the nonlocal `f` (`fn_indx`, nonlocal by design).
- **Verification** (no Professional oracle: it returns a zero field on
  the converted `mhardsoil_*` inputs): results identical for 1, 2, 3, 4
  and 8 equilibrium iterations, with and without the post-solve pass;
  hand check of the end state of `mhardsoil_plast` / `_gp0_off`:
  `sigma = sigma0 - C(E50, nu50) eps_p` exactly (zero total strain;
  HEAD violated it: sigxx -2.4885 vs -3.4962 implied by its eps_p),
  `kap = sqrt(0.5 eps_p:eps_p)`, `f = 5e-9` / `2e-7`. Tests
  `test-2026/hardsoil_kappa_iterations1/2` (RED on the old code).
  `mhardsoil_plast` (kap 0.00104991, sigxx -3.11898) and
  `mhardsoil_gp0_off` (kap 0.000494563, sigxx -1.56781) were retargeted
  as GNU-calibrated (project owner 2026-10-10).

