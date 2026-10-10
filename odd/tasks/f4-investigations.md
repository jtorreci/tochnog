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
- [x] T8 C45 + C46 (owner 2026-10-10: study, not only document): `control_mesh_copy`/`_mirror` after a time step stop with a `node_dof_start_refined` error (Pro runs); copy/switch order on the same index (Pro copies first); `control_mesh_move` has no effect in the Pro (Pro-suspect?); `mesh_activate1` final-step vely (GNU 0, Pro -9.95e-6). Route: delegated writer. Done except C46c (rotate sweep, parked: large) and two owner decisions (C46b, C46e class D).
- [x] T9 C43: GNU Bi-CG converges far slower than the Pro (`force11` >600 s vs 13 s). Route: read-only investigator + omen, owner decisions 2026-10-10 (items 1-6), then one delegated writer. Done: SuperLU UAF, interface_convert UAF, recovery table, pardiso -> SuperLU, CSR operator, docs (C46b/C46e class A, C46c not scheduled).
- [ ] T10 C48: interface tension-status record written from a double buffer (`interface_element()` ~1031) → stack garbage. Route: delegated writer.
- [ ] T11 C49: build optimisation `-O1` vs `-O2` — find why -O1 was chosen, measure -O2 (results + time) on omen; owner decision before changing the default.
- [ ] T12 C47: `force11` `mom1y_sig` 5.80223 vs Pro 5.76460 (0.65 %, model difference, solver excluded).
Order (owner 2026-10-10): T10 → T11 → T5 → T12 → T6.

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
- T8 (route: delegated writer, inline single-threaded; local probes GNU vs
  Pro 02-08-2026, scratch clean builds, omen A/B `f4t8_runs`).
  C45 `0ca144c`: mesh.cc copy/mirror lacked node_dof_start_refined
  (renumber() rc 1) and re-initialised the IP state (1 % off); now node
  items + ELEMENT_DOF duplicated, element numbers old+max+1; loaded probes
  = Pro to 4e-10. Tests mesh_copy/mirror_timestep1 RED rc 1 / GREEN.
  C46d numbering `bc2235e`: macro -bar/-rectangle/-brick x fastest (Pro
  .dbs identical); mesh_macro_numbering1/2 RED/GREEN, mesh_switch_timestep1
  check retargeted (node 3 at 0 1 = Pro node 3).
  C46a `e730052`: order copy, move, switch, mirror (pairwise Pro probes);
  mesh_order_timestep1 RED (0 2 / 4 0) / GREEN (0 7 / -1 5).
  C46b `c21b62e`: move once per index; Pro moves only node_start_refined
  (Pro-suspect, class D, owner decision pending); mesh_move_timestep1 RED
  22 1 / GREEN 12 1.
  Switch `fb834c6`: two entries = swap (manual 6.244, Pro); 3 entries GNU
  extension (class C); Pro leaves node_start_refined unswitched (C46e,
  class D, owner decision pending); mesh_switch_pair1 RED 2 0 / GREEN 0 2.
  C46d gravity `7aeac4a`: index-0-only read (input index 10), no
  bottom-to-top window, default method1, inconsistent reduced stiffness
  -> fixed; mesh_activate1 vely/veliy/sigyy = Pro to 3e-4, target met;
  mesh_activate_gravity_time1 RED rc 1 / GREEN; residual class B.
  UAF `fd477ff`: buffered node-item copies in mesh.cc (ASan, from T9);
  outputs identical; interface.cc interface_convert still has it.
  C46c parked: valid probe = Pro rotate1.dat (Pro 14 nodes/5 prism6, GNU
  rc 1: no run-time 2D->3D); rotate_angle is the sweep angle in the Pro.
  Evidence: build_safe 158 -> 160/164/166/168/170/171 OK, 0 FALLO, rc 0
  per stage (final 171). Corpus A/B (old = F4 T7 ABC = HEAD, new =
  7aeac4a): Pro other 273 -> 273 (mesh_activate1 rc 1 -> 0), sfnet 192 ->
  192, t26 99 -> 107 (+8 new), vs 232 -> 232; dbs diffs = renumbering,
  round-off / solver order (<= 5e-6), control_distribute realisations,
  pre-divergent inputs. Docs commit: see git log. Final build copied to
  build/tochnog (md5 d8efb7e1).
  Next: owner decisions C46b/C46e; C46c rotate sweep (large); interface.cc
  UAF; per-IP gravity activation (C46d residual).
- T9 (route: read-only investigation (scratch prototypes, omen) -> owner
  decisions 2026-10-10 -> delegated writer; single writer, inline edits of
  so_suplu.c, so.cc, interface.cc, miscel.cc, so_bicg.cc; clean scratch
  builds per commit, omen A/B `f4t9w_runs`, hogs `f4t9w_meas`).
  Trigger evidence: 5 code files + docs (writer trigger), corpus/hog runs
  on omen (per-action workers not needed: scripted runs).
  1 `681dc83` so_suplu.c: dQuerySpace after Destroy (UAF), L/U destroyed
  when never created, Destroy_CompCol_Matrix free() of new[] arrays.
  Test superlu_solve1 RED rc 139 (3/3) / GREEN rc 0; ASan old
  alloc-dealloc-mismatch (superlu_solve1, sfnet examp9), new clean.
  2 `789f8bf` interface.cc interface_copy_node_dbl (3 sites); ASan now
  passes interface_convert on interface_bar2_hex8 / bar3_quad8 /
  quad8_hex20; next report = new C48 (tension status written from the
  stack; also visible as run-to-run garbage in the interface3 /
  interface_quad8_hex20 .dbs).
  5 `acfac40` miscel.cc recovery factor table (quad4/quad9/hex8/hex27/
  tet10): bit-identical (297 suite .dbs, corpus old -> b1 identical apart
  from C48); large1 37 -> 27 s.
  3 `3349515` so.cc -matrix_pardiso -> SuperLU, Bi-CG fallback, RLIMIT_AS
  memory guard (70 % of the room; TN_PARDISO_SUPERLU_MB), SuperLU off
  after a failure. Tests pardiso_superlu1 / pardiso_bicg_fallback1 RED
  (old: Bi-CG, no SuperLU) / GREEN. large1 forced to pardiso under 4 GB:
  clean fallback rc 0. force11: TIMEOUT -> 20 s on omen (rc 1, mom1y_sig
  5.80223 vs Pro 5.76460: new C47, same value from the band LU).
  4 `b7009dd` so_bicg.cc CSR operator (threads, EBE fallback with mpc
  ties / TN_BICG_EBE). Test bicg_csr_threads1 RED (no CSR line) / GREEN.
  force11 ~5x iterations per second. RED/GREEN not applicable to the
  speed itself (exception: performance; checked by corpus A/B).
  6 docs commit: DIVERGENCES SOLVER-DEFAULT (B), C46b/C46e -> A, C47/C48
  (B), C46c not scheduled; BACKLOG C43 resolved, C46 decisions, new C47,
  C48, C49 (-O1 vs -O2); SEGUIMIENTO rows; solver pages; this file.
  Evidence: build_safe 171 (old) -> 175 (pardiso) -> 177 OK / 0 FALLO,
  rc 0 (main tree, final build md5 e84e554c = omen b4). Corpus A/B on
  omen, old c62d28f -> final: Pro other OK 273 -> 276 (force11,
  ground11_nonsaturated, ground11_phreatic_level finish), sfnet 192 ->
  192, t26 110 -> 111 (+superlu_solve1), vs 232 -> 232; tests finished on
  both 1180 -> 803 s (-32 %). node_dof changes: pardiso inputs <= 2.7e-6
  (Bi-CG tolerance vs direct), CSR round-off <= 1e-5 except path-
  dependent models (mesh_activate1 3.3e-3, distri3 1.8e-4, examp16
  1.4e-2, examp10 2.8e-3) where a 2-thread run of the old build moves
  them as much, or the absolute bicg_error_minimum floor binds (examp10).
  Hogs (omen, 12 GB cap, one at a time): see solver_memory.md F4 T9.
  Engram mirror: pending (writer subagent).
  Next: C47 (force11 3D moment), C48 (cheap fix + A/B), C49 (-O2),
  foundation_consolidation_settlement (Bi-CG breakdown + band retry).

