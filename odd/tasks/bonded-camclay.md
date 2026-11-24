# Bonded-soil generalised Cam clay (F5 follow-up)

Objective: implement in the GNU the three missing Professional names of the generalised non-associate Cam
clay for bonded soils (manual 6.700, 4.21, 4.38) from the published Nova model. Branch `dev/f5-keywords`.
Delivery: ask-on-risk, work-unit commits, never push. Runner: `scripts/build_safe.sh` (`TN_SKIP_SUITE=1` for
builds); single-element reference integrator in the scratchpad (not committed as a test source).

## Specs

- S1 Owner request (L1): "Intentemos implementarlo." The three names are
  `group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils`,
  `materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils_history` and
  `materi_strain_plasti_generalised_non_associate_cam_clay_for_bonded_soils`.
- S2 Model definition from the literature (Nova 2005 ICSMGE [A]; Arroyo, Castellanza & Nova 2005, S&F 45(2) [B]):
  shifted stress p* = p' + p_t, size p_c = p_t + p_s + p_m, p_t = k p_m; yield f and potential g both of the
  Lagioia-Puzrin-Potts form with their own (mu, alpha, M); hardening
  dp_s = rho_s p_s (deps_v^p + xi_s deps_s^p), dp_m = -rho_m p_m (|deps_v^p| + xi_m deps_s^p).
  Record (Q1, owner choice 1): 11 values `mu_f alpha_f M_f mu_g alpha_g M_g rho_s xi_s rho_m xi_m k`; elasticity
  from the existing GNU elasticity records; no Lode-angle dependence (one M per surface); initial p_s, p_m in the
  history dofs. Future work to document (owner): option 2 (M_c/M_e Lode dependence for f and g, 13 values) and
  option 3 (full Arroyo et al. model with Borja hyperelasticity).
- S3 History: `materi_plasti_..._history` adds history dofs `gencamhis0..` (p_s, p_m); combining it with another
  history initialisation is an input error.
- S4 Plastic strain: `materi_strain_plasti_...` adds `eppgencamxx..zz`, filled like the other per-model
  plastic strains.
- S5 Next feature after S1-S4 (L2): reopen the GiD family (62 names).
- S6 Verification: single-element tests (isotropic compression, drained and undrained triaxial, oedometer)
  against an independent integration of the same equations; the Professional cannot be the reference (it
  never ran the model, KEYWORD-GAP-REPORT 2026-11-24). Divergence row in DIVERGENCES-PROFESSIONAL.md.

## Tasks

| ID | Specs | Route | Status | Commit |
|---|---|---|---|---|
| T1 | S2 | inline | done: Q1 answered (L5) | - |
| T2 | S1, S3, S4 | inline | done: registration, history, plastic strain dofs, checks, docs (L6) | see T4 |
| T3 | S2, S6 | inline | done: plasti_rule block, closed forms, 8 tests (L6) | see T4 |
| T4 | S6 | inline | done: suite 1379 OK / 0 FALLO, probe 1076 accepted, docs | this commit |

## Log

- L1 (2026-11-24) owner: "Intentemos implementarlo. Te propondría seguir con GiD cuando hagamos estos 3. Aunque lo
  descartamos desde el principio, estoy trabajando en un proyecto de un preprocesador que admite formato GiD y nos
  puede interesar"
- L2 (2026-11-24) Professional evidence: manual gives "index ..." and one theory line (Milan model of R. Nova). Binary
  symbols gencam_set_parameters (10 doubles), gencam_fg_function, gencamhis, eppgencam*. One-element probe never ran
  (fewer than 12 values: db error; 12-20: assertion failure).
- L3 (2026-11-24) GNU map (explorer): one `plasti_rule()` in plasti.cc with a block per model; flow direction is the
  numerical gradient of f_flow; diprisco (plasti.cc:480-780) is the non-associated template. Shared history slot
  (input.cc:350-398); materi.cc pm_on[9]/pm_indx[9] must grow to 10; names (71-73 chars) fit MCHAR 100.
- L4 (2026-11-24) Literature (research worker): sources [A] https://www.issmge.org/uploads/publications/1/22/STAL9781614996569-3651.pdf
  and [B] https://www.jstage.jst.go.jp/article/sandf1995/45/2/45_181/_article/-char/en. LPP closed forms (real roots,
  complex roots, mu = 1) derived and checked numerically by the worker; Nova et al. 2003 and LPP 1996 full texts not
  accessible. [B] Table 2 gives parameter sets (LCS-1, LCL-1/2, NC-1, MCC-1) with Borja hyperelasticity and Lode
  dependence (M_c/M_e).
- L5 (2026-11-24) owner, Q1 "scope of the model": "1, y deja documentsdo el trabajo futuro para 2 y 3" (1 = core,
  11 values, GNU elasticity, one M per surface; 2 = core + Lode, 13 values; 3 = full Arroyo et al. with Borja
  hyperelasticity).
- L6 (2026-11-24) Implementation and evidence. RED: before the change the three names were rejected by the parser
  (keyword probe NO). Closed forms of I(eta) vs quadrature: 156 sets, worst rel. diff 5e-12 (harness built from the
  plasti.cc text). iso1: p = p_s + p_m and both exponential laws hold to 1e-9. tension1: yields at p = -k p_m.
  shear1: final state on the yield surface by independent quadrature (2e-9); dilatancy within ~0.01 of the model
  (first-order explicit direction). Yield-function design iterated: ln form alone fails at the apex, q form alone
  makes the critical state oscillate (spurious debonding); hybrid kept. Step study: critical ratio 1.125/1.186/1.195
  at strain steps 5e-4/1e-4/5e-5 (M_g 1.2).
- L7 (2026-11-24) Suite: 1379 OK / 0 FALLO, runs 226/226, mode000 0. Keyword probe on the new build: 1076 accepted
  (+3), 0 missing apart from the 63 discarded. Next: S5, reopen the GiD family.
