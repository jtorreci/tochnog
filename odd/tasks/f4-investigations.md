# F4 — design items and investigations

## Objective
Close the remaining design-level and unexplained differences with the
Professional, diagnosing before changing code. Expensive runs go to omen
through `scripts/remote_run.sh`.

## Tasks
- [x] T1 TPL + K10: align `groundflow_total_pressure_limit` with manual 6.587 / Pro (default 0, cap the calculated total pressure, zero the groundflow velocity where reached, never cut the `pres` dof). Owner decision 2026-10-07. Route: delegated writer.
- [x] T2 C41: `mohr_coul_direct7` (`-total_linear` + velocity `mpc_node_number` tie: unequal disy at tied nodes). Route: same writer after T1.
- [x] T3 C33: `set_stress` reads DTIME from VERSION_NORMAL (previous dt) → `_visco` direct law does not cut in the first step; principal tension cap ignores `_visco`. Route: same writer.
- [x] T4 C31 + K20: quadrature design (TET10 nodal rule degree 1; axisymmetric quad4 forced 1-point rule). Route: read-only investigator + omen A/B, then decision (owner 2026-10-09: v2 + HRZ); implementation by a delegated writer.
- [ ] T5 C35/K02: integration-point state one-iteration lag. Investigation.
- [ ] T6 Open-cause clusters: K14 (hypoplasticity history/stress path), K15 (MC/DP/von Mises), K24 (3D consolidation `large2`), K21 (Navier-Stokes), C7 (Mazars), C11 (dynamics), C3 (3D section sign), C6 (iface_gen coverage). Investigations, one by one.
- [x] T7 Memory hogs (>4 GB with any binary): Pro `force11`, `large1`, `ground19_water_under_dam`, `mesh_activate1`, `foundation_consolidation_settlement` (runs since F4 T4 K20, bad_alloc under the 4 GB cap). Route: read-only investigation, owner decision 2026-10-09 (fix A, B, C), then one delegated writer. Memory fixed; solver-speed residual of `force11`/`foundation...` open (BACKLOG C43).

## Rules
Manual/analytic decides; Pro is the oracle when the manual is silent;
divergences need an owner decision and go to DIVERGENCES-PROFESSIONAL.md.
Every tochnog run capped (`ulimit -v 4000000` + timeout, ≤4 local jobs);
corpus A/B and Pro sweeps on omen (`scripts/remote_run.sh`, ≤6 jobs).
Each fix: RED/GREEN test, build_safe green (baseline 138 OK), corpus A/B,
docs + SEGUIMIENTO + BACKLOG + DIVERGENCES.

## Progress
- Branch `dev/f4-investigations` from `af080b0`.
- T1 (route: delegated writer; inline edits in groundfl.cc, top.cc, headers):
  Pro rule implemented (default 0 via `groundflow_total_pressure_limit()`,
  total-pressure cap in `groundflow_phreatic_coord()`, velocity 0 where the
  IP total pressure reaches the limit, pres dof never cut). GNU dry
  consolidation skip removed: it conflicts with the Pro (consolidation kept
  with limit 0). Evidence: 8 gtpl probes + ground1 = Pro; build_safe
  143 OK / 0 FALLO; corpus A/B vs 231→233, t26 82→86, sfnet 172→173
  (ho_mech1 timeout flake), Pro other 242→242. Retargeted test-2014
  groundflow_total_pressure_limit/_dry (Pro values). Commit 3acce9d.
- T2 (route: same writer, inline mpc.cc): root cause mpc_node_apply
  re-synced the slave velocity after the displacement integration, so the
  slave displacement lagged one solve (-total_linear sees it). Fix:
  integrate the slave displacement with the slave velocity. Evidence:
  mohr_coul_direct7 sigxy 0.0998334 = Pro, tied nodes equal; residual
  disy/sigzz difference -> DIVERGENCES/BACKLOG C41. build_safe 144 OK /
  0 FALLO; A/B t1->t2: only mohr_coul_direct7 and the new test fixed
  (ho_mech1 timeout flake). Commit c4836d5.
- T3 (route: same writer, inline stress.cc): DTIME VERSION_NEW in
  set_stress (also the lambda*dtime of the incremental visco laws, Pro
  confirmed), visco factor on the principal tension cap, separate
  tension/MC factors in full mode (Pro probes mcstv/mcsboth). Evidence:
  14 probes = Pro; build_safe 148 OK / 0 FALLO; A/B t2->t3 only the 4 new
  tests change. No other dt mismatch in set_stress (dtime feeds only the
  visco laws and the three direct cut-offs). Commit 8a1fc75.
- T4 (route: read-only investigation with scratch variants v1-v4 and omen
  A/B, then one delegated writer; owner decision 2026-10-09: adopt v2 +
  HRZ, v4 out of scope -> BACKLOG C42).
  WU1 K20 (polynom.cc pol(), calcul_force.cc msf_element_rule): forced
  1-point axisymmetric rule and quad9 refusal removed. Tests RED old /
  GREEN new: axisym_quad4_hourglass1 (disy -4.380511327e-4 = Pro; old
  ~1e12), axisym_quad4_ring1 (disx 1.081776788e-2 = Pro, sigzz 10760.09
  vs Pro 10759.22, C42), axisym_quad9_uniform1 (1000/1000/500; old
  refused). msf_axisym moment 0 -> 1/960 (analytic). C32 closed (axisym9
  row stale, = Pro before). Commit 9d3e9f6.
  WU2 C31 (polynom.cc Keast 5-point, miscel.cc LSQ 1,L1,L2,L3 closed
  form, general.cc HRZ for principal dofs only): tet10_body_load1 sigzz
  -1.00000000515 (old -0.894), tet10_mass_hrz1 v 2.16 / 0.405 (HRZ; Pro
  -1.2 / 0.3; old 59999.99 / 0.36), tet10_dynamic_gravity1 disz -0.5286
  at t 3 (stable; row sum 5e8). Divergences TET10-HRZ (A) and
  TET10-PPOINT (B, pending owner). Commit bcbdda0.
  Evidence: clean builds; build_safe 154 OK / 0 FALLO (148 + 6), rc 0.
  Corpus A/B on omen (f4t4_runs, old 555294f -> new): Pro other 243 ->
  247 (+axisym7, force4, slide2, tet10_2; foundation_consolidation_settlement
  now runs into the 4 GB cap -> T7), sfnet 173 -> 173, t26 91 -> 94 (the
  3 new K20 tests; the 3 tet10 tests were added after the push and were
  checked by build_safe), vs 233 -> 233 (slide_axi/node_slide now a rigid
  translation, targets still met). No retargets.
- T7 (route: read-only investigation -> owner decision 2026-10-09: fix
  A, B, C -> one delegated writer; scratch clean builds, omen A/B).
  A `so.cc` dense n x n `mat_dense` (never read) removed, commit 3faf78e.
  B `so.cc` `solve_band_retry()`: the Bi-CG direct-LU retry builds its
  band matrix only after a failure, exact band, Cuthill-McKee ordering
  when narrower, commit 5bfa165. C `top.cc` `step_start_control_entry`:
  `control_mesh_macro`/`_copy`/`_mirror`/`_switch` once per control index
  (Pro probes: once), tests `mesh_macro_timestep1` (202 nodes; old 4040)
  and `mesh_switch_timestep1` RED old / GREEN new, commit ab99338.
  Measurements + docs commit 587bc13; TET10-PPOINT -> class A (owner
  decision 2026-10-09) docs commit b79afed.
  Evidence: clean builds (final = A+B+C build, md5 identical, copied to
  build/tochnog); build_safe 154 OK (A), 158 OK / 0 FALLO (A+B+C, +4 C
  checks), rc 0. Corpus A/B on omen (f4t7fix_runs, 4 GB, timeout 120):
  old -> A: only ground19 (runs) and foundation (bad_alloc -> timeout);
  A -> A+B: force11 runs (timeout), round-off in retried solves only
  (force9 4e-7, pressed_rubber_cylinder 1e-7, ho_mech1 8e-8; examp2
  bit-identical); A+B -> A+B+C: only mesh_activate1 (runs, 202 nodes,
  veliy = Pro 5e-10, final vely differs -> C46) and the 2 new tests.
  Pro other 271 -> 273, sfnet 192 -> 192, t26 99 -> 99 (+2 new), vs
  232 -> 232. Wall time (rc 0 on both) 838.6 -> 784.8 s (-6 %).
  Hogs (12 GB, one at a time, old -> new; Pro): force11 7.8 GB bad_alloc
  -> 0.27 GB, Bi-CG > 600 s (Pro 0.27 GB, 13 s); large1 bad_alloc -> rc 0
  0.65 GB 58 s (Pro 1.3 GB 9 s); ground19 5.2 -> 0.12 GB, 20 -> 18 s
  (Pro 0.17 GB 3 s, target met); mesh_activate1 11.9 GB bad_alloc ->
  0.035 GB 29 s (Pro 0.07 GB 26 s); foundation bad_alloc -> 1.5 GB,
  timeout (Pro timeout). New open items: C45 (copy/mirror after a time
  step: node_dof_start_refined error), C46 (move no-op in Pro, copy/
  switch order, rotate, mesh_activate1 vely).
