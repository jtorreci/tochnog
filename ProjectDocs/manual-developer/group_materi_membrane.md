# group_materi_membrane — developer

## Implementation

- `membrane.cc`, `membrane_apply()` (lines 24-70), called from `stress.cc`
  `set_stress()` (~line 1871) after the elastic stress (and damage), BEFORE
  the plastic-yield test.
  - 2D: `work[2] = -new_sig[8] / memmat[2][2]`, an ELASTIC zz strain
    correction that makes sigzz zero.
  - 1D: a 2x2 solve (`matrix_inverse`) for the yy and zz corrections.
  - The correction is added to `inc_ept`, `inc_epe`, `new_ept` and
    `new_epe`; `new_ept` carries the old total strain, so `eptzz` is the
    accumulated total zz strain (elastic + plastic) of the final state.
  - Returns 0 when the correction exceeds `EPS_STRAIN_CORRECTION` (1e-10);
    the caller notes it always gets 1.
- Error exits: `group_axisymmetric` in 1D, `group_materi_elasti_compressibility`
  (the check is duplicated; harmless).

## Hardcoded parameters

- `EPS_STRAIN_CORRECTION = 1e-10`.

## Verification

- `validation-suite/test-2026/membrane_eptzz_total1/2.dat` (F4 T6f,
  2026-10-12): single quad4, von Mises 24.3 with kap hardening (slope 224),
  E 7000, nu 0.2, `-total_linear`, stretched to eptxx 0.02; dt 0.01 / 0.1.
  Analytic (Tochnog kap = sqrt(epp:epp/2) = sqrt(3)/2 eppxx):
  eppxx = 0.0165286 / (1 + 224 sqrt(3)/2 / 7000) = 0.0160829,
  sigxx = 24.3 + 224 kap = 27.4199, eppyy = eppzz = -eppxx/2 (uniaxial
  Prandtl-Reuss), eptzz = -nu sigxx/E + eppzz = -0.0088249.
  Targets eptzz -0.0088249 +- 3e-4 / +- 2e-3, sigxx 27.4199 +- 0.05.
  GNU on HEAD: eptzz -0.0086300 / -0.0072102, sigxx 27.4202 / 27.4399,
  rc 0. Guard tests (no RED phase: the GNU already stores the total); the
  Professional values -9.6e-5 / -7.9e-4 fail them.
- Identity check (both dt): eptzz = -nu sigxx/E + eppzz holds to 1e-7
  (dt 0.01: -0.0007834 - 0.0078466 = -0.0086300).

## Pending refactorings

- BACKLOG C69: the correction is elastic and is applied to the trial
  stress before the plastic return, so the plastic flow of a step does not
  see an exact plane-stress state. The trace of the plastic strain stays 0,
  but its split between yy and zz is first order in dt (eppyy / eppzz
  -0.008236 / -0.007847 at dt 0.01, -0.009654 / -0.006426 at dt 0.1,
  -0.008061 / -0.008022 at dt 0.001; exact -0.0080414 each). A
  plane-stress-consistent return mapping would remove it.

## Differences with the Professional version

- **MEMBRANE-EPTZZ (deliberate, project decision 2026-10-12, F4 T6f, K15
  `validation_3`).** GNU `eptzz` = total zz strain (manual 6.685, plane
  stress; elastic + plastic). The Professional 02-08-2026 stores only the
  last step's zz correction (single element: -9.6e-5 at dt 0.01, -7.9e-4 at
  dt 0.1, analytic -0.0088). Not followed. Guard tests
  `membrane_eptzz_total1/2`. See `DIVERGENCES-PROFESSIONAL.md`.
