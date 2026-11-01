# group_materi_plasti_mohr_coul — developer

## Implementation

- `plasti.cc`, `plasti_rule()`, block `GROUP_MATERI_PLASTI_MOHR_COUL`
  (~lines 819-860, Sprint 10 `b315ade`): `matrix_eigenvalues()` on `sig`,
  sig_max / sig_min of the three principal stresses, then
  `f_yield = 0.5(sig_max-sig_min) + 0.5(sig_max+sig_min) sin(phi) - c cos(phi)`
  and `f_flow` with `phi_flow` (manual Professional 6.725).
- Flow direction: the generic central-difference gradient of `f_flow`
  (`plasti_rule()` task `GET_FLOW_RULE_GRAD`, ~lines 1375-1393; step
  `EPS_SIG * |sig|`). No special vertex/corner treatment: at a vertex the
  difference quotient averages the adjacent faces.
- The return mapping is the generic cutting-plane driver in `stress.cc`
  (`set_stress()`); `plasti_on_boundary` scales `phi` and `phi_flow`.

## Hardcoded parameters

- `EPS_SIG` (finite-difference step of the flow gradient, relative to the
  stress norm).

## Verification

- `validation-suite/test-2014/mmc_tension.dat`: sigxx 46.19 =
  2c cos(phi)/(1+sin(phi)) (analytic).
- `validation-suite/test-2026/mmc_tension_vertex1/2.dat` (F4 T6f,
  2026-10-12): the same model at dt 0.05 / 0.047 with targets sigxx 46.19
  +- 0.2, sigzz 0 +- 0.05, eppzz 0 +- 1e-4. GNU on HEAD: sigzz -3.3e-10 /
  1.1e-9, eppzz 3.3e-13 / -1.1e-12 (rc 0). Guard tests: GNU already
  matched the continuum answer, so no RED phase applies. The Professional
  02-08-2026 values (sigzz 0.974, eppzz -9.7e-4 at dt 0.05; sigzz 0.367 at
  dt 0.047) fail them.

## Pending refactorings

- None specific to this law. A consistent corner/vertex return (Koiter
  multi-surface) would make the vertex exact by construction instead of
  relying on the central-difference gradient.

## Differences with the Professional version

- **MC-VERTEX-ZZ (deliberate, project decision 2026-10-12, F4 T6f, K15).**
  `mmc_tension`: GNU sigzz ~1e-9 / eppzz ~1e-12 for every dt (continuum
  solution at the vertex sigyy = sigzz = 0). The Professional 02-08-2026
  gives sigzz 0 / 0.974 / 0.367 / 1.44 for dt 0.1 / 0.05 / 0.047 / 0.0047
  (eppzz -9.7e-4 at dt 0.05): an algorithmic vertex artefact, not a
  property of the law. Not followed. Guard tests `mmc_tension_vertex1/2`.
  See `DIVERGENCES-PROFESSIONAL.md`.
