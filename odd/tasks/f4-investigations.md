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
- [x] T5 C35/K02: integration-point state one-iteration lag. Investigation, owner decision 2026-10-10 (post-solve material state pass, gate inertia/slide/contact; mhardsoil blocker: option 1, hardsoil kappa fix first). Done `8e4d080` + `e8e3f12`.
- [ ] T6 Open-cause clusters: K14 (hypoplasticity history/stress path), K15 (MC/DP/von Mises), K24 (3D consolidation `large2`), K21 (Navier-Stokes), C7 (Mazars), C11 (dynamics), C3 (3D section sign), C6 (iface_gen coverage). Investigations, one by one.
- [x] T7 Memory hogs (>4 GB with any binary): Pro `force11`, `large1`, `ground19_water_under_dam`, `mesh_activate1`, `foundation_consolidation_settlement` (runs since F4 T4 K20, bad_alloc under the 4 GB cap). Route: read-only investigation, owner decision 2026-10-09 (fix A, B, C), then one delegated writer. Memory fixed; solver-speed residual of `force11`/`foundation...` open (BACKLOG C43).
- [x] T8 C45 + C46 (owner 2026-10-10: study, not only document): `control_mesh_copy`/`_mirror` after a time step stop with a `node_dof_start_refined` error (Pro runs); copy/switch order on the same index (Pro copies first); `control_mesh_move` has no effect in the Pro (Pro-suspect?); `mesh_activate1` final-step vely (GNU 0, Pro -9.95e-6). Route: delegated writer. Done except C46c (rotate sweep, parked: large) and two owner decisions (C46b, C46e class D).
- [x] T9 C43: GNU Bi-CG converges far slower than the Pro (`force11` >600 s vs 13 s). Route: read-only investigator + omen, owner decisions 2026-10-10 (items 1-6), then one delegated writer. Done: SuperLU UAF, interface_convert UAF, recovery table, pardiso -> SuperLU, CSR operator, docs (C46b/C46e class A, C46c not scheduled).
- [x] T10 C48: interface tension-status record written from a double buffer (`interface_element()` ~1031) → stack garbage. Route: delegated writer. Done `c06fa2d`.
- [x] T11 C49: build optimisation `-O1` vs `-O2` — find why -O1 was chosen, measure -O2 (results + time) on omen; owner decision before changing the default. Study done (default unchanged); owner decision pending.
- [x] T12 C47: `force11` `mom1y_sig` 5.80223 vs Pro 5.76460 (0.65 %, model difference, solver excluded). Root cause: quad9 -> hex27 extrusion did not share the inner-plane nodes. Done `24301c8`.
- [x] T13 C49 prerequisite: audit every `db()` / array call that receives a fixed-size scratch buffer (`idum[1]`, `ddum[1]`, small arrays) and can write more values; fix them; run build_safe + corpus under ASan; then switch the build to `-O2` with a full corpus A/B (owner 2026-10-10). Done `40ae7cd` (fixes + audit script) + `44c68c7` (`-O2` default).
Order (owner 2026-10-10): T10 → T11 → T5 → T13 → T12 → T6.
T6 triage (2026-10-11, HEAD 53dbb86) and owner order (2026-10-11):
- [x] T6a Cheap batch: C55 (`print_un.cc:149` node_rhside index), C56 (`input.cc` `cin >> str` overflow), `-size_tot` (calcul.cc sqrt(0.5 a:a) vs Pro sqrt(a:a); manual first), C3 (thickness_switch default; sign = plot_switch convention), K14a (`hypoplas.cc` hyhis write-back), C6 (corpus input repair). Done `f7559bc` `80b19af` `4705692` `7e09c40` `523a959` `3e7837b` (+ docs); C3 square-face tie -> owner decision.
- [x] T6b K24 `large2` transient 3D consolidation coupling (4.5 % pres, dt-independent; drained steady state equal). Closed 2026-10-11: deliberate divergence class A (`c350b3c`, `0cc417f`).
- [x] T6c K15 `validation_14` `control_restart` semantics (without restart every node_dof = Pro). Done `8ca6124`: restart resets the IP state (manual 6.368, Pro deletes element_dof); `validation_14` NOT caused by the restart (re-scoped, open under T6f); element history records -> BACKLOG C58.
- [x] T6d C11 `dynamic3` (disp_primary_gate per material) -> C52. Done `28992d5` (C11), `e0f2204` (C59), `b318484` (interface history), `c276931` (inertia pass).
- [x] T6e C53 slide in the post-solve pass. Done `b4e4fa4` (slide law), `07433f5` (slide/contact in the pass, gate removed).
- [ ] T6f Rest by relevance (K15 `validation_14` done `1d446b3`; K15 validation_3 kap done `c3e6d2a`; C68 IP states kapsh/cap1/damage/rho/Maxwell done `70235c6`; C68d Maxwell old-state index done `2327083`; K21 recovery part done `2305659`, gate removed `e2b044f`; `-ra` ranges done `35025dc`): K15 `validation_14` (re-scoped by T6c: pass-1 hole-boundary nodal stress/plastic strain up to ~20 % at 8 nodes on the r=0.01 circle -> `refine_locally` picks a different element -> pass 2 on a different mesh), K15 mmc_tension and validation_3 eptzz closed 2026-10-12 as class A divergences (`600a0b4`), K14b hypo3 and C7/C7b Mazars closed 2026-10-12 as class A divergences (`ae0691b`), K21 (recovery done; validation_1 stokes default + convection remain), C57 (hex27 nodal stress resolved by C70; post_point open), C70 done `39d0138` (all Gauss rules + 3D physical fit), C71 (tube1 -lobatto NaN), C75 (prism6 points). Parked: K14c (ISA, class D), mohrcou2 sfnet (input artifact).
- [x] T6f C72 convected nodal stress + validation_1 steady guard (owner-approved 2026-10-13; route: delegated writer, trigger: general.cc + 3 tests + 9 docs). Done `c1a1ae9` (fix + tests), `cc245f7` (docs). Tests registered in `scripts/build_safe.sh` `421a3f4`.
- [x] T6f C70 nodal recovery of the other Gauss rules and the 3D elements in physical coordinates (route: delegated writer; trigger: multi-file code + 7 tests + build_safe + 9 docs). Done `39d0138` (fix + tests + registration), `170e573` (docs). Found C75 (prism6 integration points).
- [x] T6f C75 prism6 integration rule: investigation + owner decision + docs (route: delegated writer; trigger: Pro/GNU probes + 7 docs). Gate stop before code; owner 2026-10-13: keep the GNU interior rule (class A). Done `4ab9549` (docs). Found C76 (`-minimal` ignored for prism6).

## Later (owner 2026-10-10, after convergence with ALL divergences documented and accepted)
- F6: exhaustive profiling and optimisation pass; special attention to the hot
  helpers used everywhere (`db()`, `pri()`): leaks, per-call costs (e.g. `pri()`
  flushes `cout` on every call), buffer contracts.

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
- T10 (route: delegated writer, inline single file interface.cc + test +
  build_safe registration + docs; clean scratch build from `git archive`,
  Pro probes local (4 GB, <=4 jobs), corpus A/B on omen `f4t10_runs`).
  Trigger evidence: one non-trivial code file (inline), corpus on omen.
  Root cause: `rec_status` was a `double*` over `get_new_int()` and the
  PUT passed `idum[1]` as the integer array (length n_intpnt): stack
  garbage (`2 -total_linear 9 0 1406...`, changing between runs). Fix
  `c06fa2d`: `long int` array of `-CLOSED`/`-OPENED`, PUT as integers,
  `delete[]` for the five output buffers. Nobody reads the record back
  (`print_only`; top.cc only sizes/allocates it).
  RED/GREEN: new test-2026 `interface_tension_status1` (two quad6
  interfaces, one closed, one open behind a gap): Pro `-closed` x3 /
  `-opened` x3; old GNU `2 -total_linear 3` (RED); new = Pro, two runs
  byte-identical (GREEN). Pro sweep, 20 `interface*` inputs of Pro
  `test/other`: all Pro records `-closed`, GNU line-equal except
  `interface_bar2_hex8` (element numbering after extrude+convert: new
  BACKLOG C50 / DIVERGENCES row). Determinism: two GNU runs identical
  on all 20 apart from `time_calculation`.
  Evidence: build_safe TN_CHECK_ONLY 181 OK / 0 FALLO, rc 0 (177 + 4
  new checks); build md5 25040c72 (= build/tochnog = omen f4t10_new).
  Corpus A/B old e84e554c -> new: classes unchanged (Pro other 276,
  sfnet 192, t26 112, vs 232); only this record changes (28 + 25 + 4
  files), plus 2-thread round-off in bicg_csr_threads1 (6.7e-11).
  Docs: user/dev group_interface, SEGUIMIENTO (row + checklist),
  BACKLOG C48 resolved + C50, DIVERGENCES C48 row removed / C50 added.
  Engram mirror: pending (writer subagent).
- T11 (route: delegated writer, study only; scratch clean builds with
  `TN_EXTRA_FLAGS` (no repo script edit), corpus A/B and timing on omen
  `f4t11_runs` / `f4t11_meas`). Builds from HEAD + T10 (O1 md5 25040c72,
  O2 08720f95, O2 -fno-strict-aliasing 6316eb79, O3 1001a41c).
  Why -O1: no recorded reason (2026-08-02 modernisation command copied
  into build_safe.sh 73f0e70; upstream Linux targets use -O2).
  Build: 35/44/45/49 s, largest compiler process 137/149/151/158 MB.
  Times (median of 3): large1 23.4/16.9/22.2/16.2 s, force11
  14.5/8.4/14.0/7.7 s. Corpus: all finished .dbs byte-identical to O1;
  new RUNFAIL at O2/O2nsa/O3: Pro interface7, interface10, mesh_cut_2,
  vs creset_iface; build_safe 177 OK / 4 FALLO (bmom* cut distance 0,
  creset_iface). Cause: latent 1-element scratch-buffer overflows (C48
  class) in print_beam_force_moment.cc, print_interface_stress.cc,
  delete.cc, top.cc. Recommendation: keep -O1; fix the buffer class
  (ASan), then -O2 with a full A/B. Details: BACKLOG C49.
- T5 (route: delegated writer; inline top.cc/elem.cc/stress.cc/headers +
  5 tests + build_safe registration + docs; clean scratch builds from
  `git archive` / `git stash create`; probes local (4 GB, <=4 jobs),
  corpus A/B on omen `f4t5w_runs`). Trigger evidence: few code files,
  investigation already mapped (scratch prototype TN_C35).
  Owner decision 2026-10-10: post-solve material state pass as the
  default (no env), no matrix store, gated OFF for inertia / slide /
  contact (temporary divergence class B `POSTPASS-GATE`, dynamic3 /
  slide4; BACKLOG C52/C53), widening only after the pass covers that
  physics + full A/B. Blocker found: test-2014 mhardsoil_plast /
  _gp0_off moved (kap 0.0006 -> 0.0030) because the hardsoil yield read
  the nodal kappa of the previous evaluation (period-2 oscillation in
  the iteration count, already in HEAD). Owner 2026-10-10 option 1: fix
  first.
  Work unit 1 `8e4d080` fix(plasti) C51: set_stress passes pl_unknowns
  (kappa / kappa_shear = step-start + own return increment) to the yield
  rules (hardsoil, kappa-hardened MC, von Mises Nadai). Sweep 1/2/3/4/8
  iterations, with and without the pass: identical (old 0.0030 odd /
  0.0006 even). Hand check: sigma = sigma0 - C(E50,nu50) eps_p exact,
  kap = sqrt(0.5 eps_p:eps_p), f 5e-9. Tests hardsoil_kappa_iterations1/2
  RED (HEAD rc 1) / GREEN. Retarget (owner, GNU-calibrated, no Pro
  oracle): mhardsoil_plast kap 0.00104991 / sigxx -3.11898,
  mhardsoil_gp0_off kap 0.000494563 / sigxx -1.56781. build_safe (commit
  1 tree) 183 OK / 0 FALLO, runs 226/226. Corpus head -> hsA: classes
  unchanged, farther 0 (changes: tensile_instability, taylor3, examp19,
  examp3, ho_mech2, examp14, mhardsoil x2; no Pro value or ~same). Same
  pattern not changed (cap epp_vol, cap1 pc, void fraction): BACKLOG C51.
  Work unit 2 `e8e3f12` feat(top): top.cc pass + post_solve_state_pass_gated,
  elem.cc matrix-store skip, headers. Probes all = Pro; tests
  post_solve_state_stretch1/_plastic1/_rigid1 RED (HEAD rc 1) / GREEN.
  Matrix vs no-matrix: probes and 742 non-gated corpus tests identical
  (bicg_csr_threads1 2-thread round-off aside). build_safe 187 OK /
  0 FALLO, runs 226/226, rc 0 (final md5 2c825fe4 = build/tochnog = omen
  f4t5w_fin2; HEAD 25040c72). Corpus HEAD -> final: closer 11, farther
  0, 94 gated tests unchanged (dynamic3, slide4), classes unchanged
  except ground11_nonsaturated OK -> TIMEOUT (84 -> 105 s alone; Pro
  times out); wall +8.8 / +7.2 / +4.7 / +12.0 % (Pro other / sfnet / t26
  / vs). Docs: manual-*/post_solve_material_state_pass.md (+ README),
  nodal_stress_post_solve_sync notes, SEGUIMIENTO, BACKLOG C35 resolved
  + C52/C53, DIVERGENCES POSTPASS-GATE (C35 row removed), AUDIT K02.
  Engram mirror: see report.
- T13 (route: delegated writer; trigger evidence: 21 source files touched,
  mechanical per-site fixes found by a script, so one writer; clean scratch
  builds from `git archive` / `git stash create`, local ASan build_safe,
  corpus A/B + ASan corpus + timing on omen `f4t13_runs` / `f4t13_meas`).
  Audit: `scripts/audit_db_buffers.py` (db() reads vs `data_length`,
  get_group_data, constant-length array_*, geometry()/msf_element_rule()
  contracts). HEAD dca225a: 3450 calls, 45 hits at 34 call sites (26 db/
  array, 7 geometry() callers with ddum[1], msf nper[2]); all fixed; 0
  hits after (exit 0). False positives met while building the script
  (prototype `sig[6]` matched as a global, unbounded maxwell length) were
  fixed in the script, none left. 38 run-time-record calls reviewed by
  hand -> C54. Work unit 1 `40ae7cd` fix(db).
  ASan+UBSan (-g -O1): build_safe head 28 failing runs / 18 FALLO ->
  fix 1 / 1 (dofrhside, C55); corpus head 45 ASan tests at 9 sites ->
  fix 2 (dofrhside C55, Pro ground19_water_under_dam parser token C56);
  no UBSan report. -O1 fix (md5 db9254c8) vs HEAD (2c825fe4): corpus
  .dbs identical (bicg_csr_threads1 2-thread round-off), classes equal;
  build_safe 187 OK / 0 FALLO, 226/226.
  -O2 (fix2 md5 cc867346 = final default build): all finished .dbs
  byte-identical to -O1; creset_iface, bmom*, interface7, interface10,
  mesh_cut_2 run; alternate1, cylinder_buckle, ground11_nonsaturated
  TIMEOUT -> OK (Pro other 275 -> 278 OK); build_safe 187 OK / 0 FALLO
  rc 0. Time (omen, alone, median of 3) large1 31.4 -> 22.2 s, force11
  18.8 -> 9.9 s; corpus wall (finished on both) Pro other 696.5 -> 582.5
  s, sfnet 240.1 -> 220.1, t26 35.9 -> 32.1, vs 28.2 -> 26.0; peak RSS
  same. Work unit 2: build_safe.sh -O2 (+ makefile note, QUALITY-CI,
  AGENTS, BACKLOG C49 resolved + C54-C56, SEGUIMIENTO row, dev page
  db_buffer_contract) = `44c68c7`. Final clean build from 44c68c7
  (default -O2): md5 019e3620 = build/tochnog (differs from cc867346 only
  in the __DATE__ string, 21 bytes); build_safe TN_CHECK_ONLY 187 OK /
  0 FALLO, runs 226/226, rc 0. Not changed: scripts/install_superlu.sh
  still builds SuperLU with -O1 (outside the edit surface)..
- T12 (route: delegated writer, single source file + test + docs; root
  cause from a prior read-only investigation, scratch `f4t12/`). Fix
  `mesh.cc` `mesh_extrude()` quad9 -> hex27: plane q (1..2nq-1) of source
  node inod = `next_node + (q-1)*nbase + inod`, created once, shared
  (numbering kept GNU-style, not the Pro compact one: renumbering would
  move the boundary-copy block of the linear branches). Clean scratch
  builds (`git archive` a08992b md5 019e3620 = old; `git stash create`
  tree md5 edb81804 = new = build/tochnog, sources = 24301c8).
  force11 (Pro 02-08-2026, `-matrix_iterative_bicg` for the Pro run):
  nodes 15600 -> 12000 (Pro 12000), rc 1 -> 0, mom1y_sig 5.80223 ->
  5.764609 (Pro 5.764612), nodal displacements = Pro to 5.5e-7, section
  forces 2e-4; nodal stresses up to 3.38 % of the field maximum off
  (largest at the loaded points x 0 y +-50.1: sigxx -816.80 vs -789.09;
  ratio 0.987..1.035, mean 1.0001; sigxy 0.9 %) -> new BACKLOG C57 (+
  `post_point` on a hex27 node not equal to the nodal value: disy
  -3.51050e-5 vs -3.51802e-5), DIVERGENCES C57 (class B).
  Test `extrude_quad9_shared_nodes1` (2 quad9 x 2 hex27 layers,
  cantilever plate; targets = Pro nodal displacements at the same
  coordinates): RED old rc 1, 84 nodes (planes 15/18/18/18/15); GREEN rc
  0, 75 nodes (15 per plane) = Pro, displacements = Pro to 2e-12.
  build_safe TN_CHECK_ONLY: 189 OK / 0 FALLO, runs 226/226, rc 0
  (baseline 187 + run_ok + node-count check). audit_db_buffers.py exit 0.
  Corpus A/B on omen (`f4t12w_runs`, 6 jobs, 4 GB, timeout 120): every
  control_mesh_extrude input (Pro other extrude1/2, safety2,
  interface_bar2_hex8, force11, Pro validation_14_mesh, new test): only
  force11 (rc 1 -> 0) and the new test change; the others keep their
  class (extrude1/safety2/validation_14_mesh PARSE, extrude2 RUNFAIL,
  interface_bar2_hex8 OK) with identical .dbs. Full Pro other (330):
  .dbs differ only for force11; classes equal except cylinder_buckle
  OK 119.0 s -> TIMEOUT 120 s (borderline; rerun alone 80 / 81 s, .dbs
  identical but time records). Work unit `24301c8` fix(mesh) (fix + test
  + build_safe + docs: control_mesh_extrude user/dev, SEGUIMIENTO row,
  BACKLOG C47 resolved + C57, DIVERGENCES C47 closed + C57, solver.md /
  matrix_pardiso.md notes). Next: T6.
- T6a (route: delegated writer, single writer; inline edits of print_un.cc,
  input.cc, calcul.cc, hypoplas.cc, corpus inputs; local Pro probes, clean
  scratch build md5 d1aa27b8 = build/tochnog = omen f4t6a_new; omen A/B
  `f4t6a_runs`). Trigger evidence: 4 code files + docs (writer trigger).
  C55 `f7559bc`: `node_rhside[ipuknwn]`; Pro probe with derivatives =
  new GNU (vely_rhside +-1.0416667, sigyy_rhside -2.0833333; old 0 and
  4.7e-310). Test dof_rhside_derivatives1 + file check (RED/GREEN).
  C56 `80b19af`: `input_word()` (setw(MCHAR)) at all 16 word reads;
  split words allowed in comments, `Word too long` elsewhere. Tests
  input_long_comment_word1 (rc 0), input_long_word_error1 (old SIGSEGV).
  range.cc still has 2 `cin >> str` (outside the surface).
  size_tot `4705692`: manual post_calcul = sqrt(Aij Aij) / sqrt(Bij Bij),
  no 0.5 (also -size_dev); Pro calcul1 1 / 0.816497. Test
  post_calcul_size_tot1 (RED/GREEN). A/B: Pro ctvmis2 1 -> 0; sfnet
  calcul1/2, ctvmis2, examp18 0 -> 1 (stale old-convention targets,
  BACKLOG B14).
  K14a `7e09c40`: hyhis1/2/3/5/7 written back (clay from statev + L;
  basic phi_mob, OCR p_e/p without s, sqrt(M:M) of the zero-increment
  tangent). hypo7/8/12/13 hyhis2/3/5 = Pro to ~3e-4; residuals hyhis1
  (GNU substep vs Pro dt), hypo12 hyhis3 (m_R L), hypo9 -> DIVERGENCES
  K14. Test hypo_masin_history_slots1 (RED/GREEN). A/B: only hyhis
  columns change. hypo7 1.7 -> 3.1 s user.
  C6 `523a959`: generated interface elements removed from the 4
  iface_gen inputs (+ element_geometry tags in iface_gen_method); runs
  generate them; 4 build_safe checks.
  C3 `3e7837b`: no code change; gforce10 input regains the force10
  plot_switch/thickness_switch -yes, target -5000 = Pro exactly;
  thickness_switch default -no in manual, Pro and GNU. OWNER DECISION:
  square end face without -yes (Pro degenerate frame vs GNU tie by the
  reference point).
  Evidence: build_safe 199 OK / 0 FALLO, runs 226/226, rc 0 (baseline
  189); audit_db_buffers exit 0. Corpus A/B old edb81804 -> new: Pro
  other OK 277 -> 277 (ctvmis2 exit 1 -> 0; dbs diffs calcul1/2,
  ctvmis2, hypo7/8/9/12/13), sfnet 192 -> 192 (exit 0 -> 1 calcul1/2,
  ctvmis2, examp18; size x sqrt(2) in examp3/10/14/17/23), t26 122 ->
  122 (new tests GREEN; bicg_csr_threads1 2-thread round-off), vs 232 ->
  232 (dbs identical). Docs commit: see git log.
- T6a C3 owner decision 2026-10-11: keep the GNU square-face tie-break, class A (`a0cea2e`). T6a closed. Next: T6b K24.
- T6b K24 (2026-10-11), route: delegated writer (docs in 11 files + 4 test
  inputs + build_safe wiring; preparation reading delegated). Cause:
  `solver_matrix_symmetric -yes` - the Pro symmetrizes the non-symmetric
  Biot u-p matrix (manual 6.1052), the GNU keeps it (Bi-CG, warning).
  Owner decision option A: GNU kept, DIVERGENCES K24 -> class A (A 20,
  B 38), also Pro-suspect. Guard tests `c350b3c`
  (test-2026/k24_consol_3d1/2 -yes/-no, k24_consol_3d_undrained1/2),
  targets = Pro 02-08-2026 `-no` re-run here at node 1: pres -0.7085296
  / disx 0.05420111 (GNU -0.708335 / 0.0542098) and pres -0.847677 /
  disx 0.01748967 (GNU -0.8476521 / 0.0175445); GNU identical for -yes
  and -no; Pro -yes -0.678332 / -0.803176 fails them; a wrong target
  gives rc 1. Guard tests: GREEN on the current binary, no RED (no
  behaviour change). Docs `0cc417f`: DIVERGENCES, AUDIT K24 closed,
  SEGUIMIENTO Descarte + register, solver_matrix_symmetric user/dev
  (manual number corrected to 6.1052), solver.md stale "-yes forces CG"
  text fixed, groundflow_consolidation_apply notes. Evidence: build_safe
  TN_CHECK_ONLY=1 203 OK / 0 FALLO, runs 226/226. The two stale
  so_bicg.cc comments were fixed afterwards (parent, comment-only).
  Next: T6c K15.
- T6c K15 `control_restart` (2026-10-11), route: delegated writer (one
  code file restart.cc + 2 tests + build_safe + docs; probes local,
  corpus A/B on omen `f4t6c_corpus`). Cause: `restart()` reset node,
  node_dof and time_current but left `ELEMENT_DOF_INITIALISED` = 1, so
  the old IP state survived (manual 6.368: restart from the start; Pro
  02-08-2026 deletes element_dof/element_intpnt_dof, sets
  number_timesteps 0). Fix `restart_reset_element_state()` clears the
  flag (NORMAL/NEW). GNU has no number_timesteps record (noted). Tests
  test-2026/control_restart_reset_ip1 (DP) / _ip2 (elastic), Pro targets
  at node 4: RED old binary rc 1 (sigyy +19.65 / +52.08), GREEN new
  (sigyy -24.18334 / -52.08333). Probes vs Pro: truss force kept by
  both (1.5), undrained pressure = Pro (-9.75297); interface, spring,
  contact-spring, tendon, damage effective stress not verified -> BACKLOG
  C58 + DIVERGENCES C58 (class B 39, total 76). Corpus A/B old d1aa27b8
  -> new 929f446a (omen, Pro reference): every control_restart input
  (`examp6`, `examp23`, Pro `validation_14`, Pro other
  `slope_nonlocal_refine`) .dbs identical except time_calculation; exit
  codes unchanged (examp6 rc 1 both = stale target; slope_nonlocal_refine
  timeout 300 s in old/new/Pro); new tests rc 1 -> 0. K15 `validation_14`
  re-scoped (AUDIT, DIVERGENCES K15), kept open under T6f. Evidence:
  build_safe TN_CHECK_ONLY=1 205 OK / 0 FALLO, runs 226/226;
  audit_db_buffers exit 0. Work unit `8ca6124`. Next: T6d C11.
- T6d C11 (2026-10-11), route: delegated writer (code in general.cc,
  so.cc, truss.cc, materi.cc, elem.cc, headers; 6 tests + build_safe +
  docs; probes local, corpus A/B on omen `f4t6d_runs`). WU1 C11 `28992d5`:
  `materi_dynamic_effective()` (control_materi_dynamic at ICONTROL
  overrules materi_dynamic, manual 6.141) in the displacement-primary gate
  and every reader, no PUT/restore of MATERI_DYNAMIC; truss default memory
  `-total_linear` with / `-updated_linear` without materi_displacement
  (Pro probes; manual 6.777 silent); signed projected length in the linear
  branch. Tests test-2026 truss_memory_default_displacement1,
  truss_total_linear_inversion1, control_materi_dynamic_only1/_half1,
  truss_dynamic_no_record1, truss_memory_default_velocity1: RED old rc 1,
  GREEN = Pro to 1e-10. build_safe TN_CHECK_ONLY=1 211 OK / 0 FALLO, runs
  226/226; audit_db_buffers exit 0. A/B old 929f446a -> c11 217e574c: Pro
  other 276 -> 277 (dynamic3 = Pro; alternate1 timeout -> rc 1), sfnet
  192 -> 192 but truss2/truss3 exit 0 -> 1 (stale targets of the legacy
  `-updated` default, Pro = new GNU -> BACKLOG B15), test-2026 only the
  new tests, test-2014 232 -> 232 (bmom_truss 2 %, targets met). Corpus
  check of the Pro memory input errors: sfnet truss5 hit -> owner decision
  2026-10-11: replicate both checks (C59, next work unit), then C52.
  WU2 C59 `e0f2204` (owner decision 2026-10-11): truss.cc input errors
  for explicit `-updated`/`-updated_linear` + materi_displacement and
  `-total`/`-total_linear` without it (GNU message names both records;
  check_error -no skips). Corpus scan: only sfnet truss5 hit (a) ->
  adapted (external corpus + omen copy; record dropped, target -93.8 =
  Pro run of the adapted input; validation-suite/test.zip copy untouched).
  Tests truss_memory_updated_displacement_error1 /
  truss_memory_total_velocity_error1: RED old rc 0, GREEN rc 1 + message.
  build_safe 213 OK / 0 FALLO, runs 226/226; audit exit 0. A/B c11 ->
  c59 20b5dc20: only the two new tests (expected error) and the
  alternate1 timeout flake; no .dbs change.
  WU3 C52 (inertia part of `post_solve_state_pass_gated`): built c52
  b4f0f5eb on top of C59; A/B c59 -> c52 NOT clean (test-2014 iface_gen,
  iface_gen_geom, iface_gen_method exit 0 -> 1; build_safe 223/226 runs)
  -> top.cc not committed, gate kept, BACKLOG C52 updated. T6d stays open
  for C52 (interface + inertia in the pass). Next: T6e C53 or C52 rework.
  WU4 C52 interface (2026-10-11, owner decision), route: delegated writer
  (interface.cc + comments in data.cc / print_interface_stress.cc, tests,
  docs). Root cause of the iface_gen* break: the interface MC history
  stored the elastic trial kt*gamma_total (no unloading after a slip
  reversal). Fix: history = clamped force. Test-2026
  interface_mc_slip_reversal1: RED old rc 1 (5.0), GREEN 4.0 (Pro
  4.000000083). test-2014 iface_mc_slip/_1step/iface_mc_num original 2014
  targets 0/0/50 restored + intpnt shear target (old rc 1). data.cc
  probe: strain reset Pro 4.0 = new (old 5.0), full reset -1.0 all;
  no code change. build_safe 214 OK / 0 FALLO, runs 226/226; audit
  exit 0. A/B c59 -> wu1 ea7181d3: only the new test (OK) and the
  alternate1 timeout flake; Pro interface*/mohr_coul_direct <= 8e-12.
  Found C60 (intpnt normal sign vs node order). Commit `b318484`.
  WU5 C52 gate (route: same writer, top.cc + docs): inertia removed from
  post_solve_state_pass_gated(). RED c52 iface_gen/_geom/_method rc 1,
  GREEN rc 0. build_safe 214 OK / 0 FALLO, runs 226/226; audit exit 0.
  A/B wu1 -> wu2 f6cb4e19: classes/exit codes unchanged (alternate1
  flake back to timeout); .dbs only in inertia models (taylor3 4.1e-4,
  distance to Pro unchanged; iface_gen* <= 0.15, targets met). ifx
  probe (interface + inertia) gap vs Pro independent of the pass ->
  C61. Commit `c276931`. T6d closed. Next: T6e C53.
  B15 (owner-approved 2026-10-11, route: delegated writer, corpus +
  docs only, no source change): sfnet truss2/truss3 adapted as the
  Professional 02-08-2026 copies other/truss2, other/truss3: explicit
  `group_truss_memory 0 -updated`; truss3 also replaces
  materi_displacement by materi_velocity_integrated (as Pro; avoids the
  C59 error). External copies external-downloads/sfnet/extracted/test/ +
  omen f4t4_corpus/sfnet (validation-suite/test.zip copies untouched).
  Original targets kept: GNU rc 0 both, node 1 rhside (0.100001,
  0.0999992); Pro own copies and the dialect-translated adapted inputs
  (0.1000002, 0.0999998), rc 0. Docs: BACKLOG B15 left for the parent
  (outside the writer surface).
  T6e C53 (2026-10-11, owner-approved), route: delegated writer (trigger:
  2+ non-trivial source files slide.cc + top.cc, tests, docs; patch from
  the read-only investigation, A/B on omen f4t6e_runs).
  WU1 `b4e4fa4` fix(slide): sig_t default 1.e300 (absent
  slide_plasti_tension = no tension limit), cap c + mu*|S| (Pro),
  slide_penalty_law principal-dof index vel_indx/nder+idim. test-2014
  node_slide recalibrated (owner): node_rhside 4 vely -5e13 (analytic
  -1e15*(v.n)*dt*n_y, rel 1e-6) + node_dof 4 eptxx 0; RED HEAD rc 1
  (vely 2.97e-12, eptxx -10), GREEN rc 0 (CORPUS-STALE-TARGETS 7).
  build_safe TN_CHECK_ONLY=1 214 OK / 0 FALLO, runs 226/226.
  WU2 `07433f5` feat(top): pass calls parallel_contact + slide();
  post_solve_state_pass_gated() deleted (both headers); DIVERGENCES
  POSTPASS-GATE removed (B 40 -> 39, total 76). Tests test-2026
  slide_friction_cohesion1 (RED HEAD rc 1 -0.002, slide fix alone rc 1
  -2.066157e-3, GREEN -2.06548203e-3 vs Pro -2.065482222595e-3) and
  contact_block_pass1 (RED HEAD rc 1, GREEN; Pro node_dof <= 7e-13).
  Pro slide4 rc 0 = Pro (HEAD rc 1). Focus set slide/contact 18/18 rc 0.
  build_safe 216 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0.
  BACKLOG C53 resolved; added C62 (pass skips parallel_new_dof_before),
  C63 (slide4 post_node_rhside_fixed vs Pro, pre-existing). No full
  corpus A/B re-run by the writer (investigation A/B covered the same
  patch). Next: T6f.
  T6f WU1 K15 `validation_14` (2026-10-11, owner-approved), route:
  delegated writer (trigger: source + 6 tests + docs; probes and A/B on
  omen `f4t6f_w`, `f4t6f_runs`, `f4t6f_val`). `1d446b3`
  fix(recovery): tria3 (1 pt and 3-pt nodal), tria6, tet4 1 pt, bar2
  1 pt nodal stress/strain recovery = arithmetic mean (Professional,
  two-element probes exact; GNU area-weighted). RED/GREEN: test-2026
  tria3_nodal_recovery_mean1/2, tria6/bar2/tet4_nodal_recovery_mean1,
  test-2014 validation_14_pass1 (rc 1 HEAD, rc 0 fix). Full Pro
  validation_14 = Pro within 1.2e-4 (HEAD velx 99 %). Corpus A/B head
  vs wu1: rc changes only the 6 new tests; closer to Pro 2 (axisym3,
  force17), farther 0; sfnet examp23 different mesh (no Pro, rc 0).
  build_safe TN_CHECK_ONLY=1 222 OK / 0 FALLO, runs 226/226;
  audit_db_buffers exit 0. BACKLOG C64 (quad4 operator), C65 (Pro
  post_point at a node), C66 (tria6 element 1e-5), C67 (unmeasured
  rules). Next: WU2 (integration-point kappa).
  T6f WU2 integration-point kappa (2026-10-11), route: delegated writer.
  STOPPED before commit (rule: a corpus target encoding the old
  behaviour -> owner decision). Fix in the working tree (uncommitted;
  patch copy in the session scratchpad t6fw/wu2_pending.patch): elem.cc
  ELEMENT_DOF read/write masks + kap_indx, materi.cc kappa block outside
  the h guard with strain_recovery_weight + IP write, general.cc kappa
  denominator weight. RED/GREEN: test-2026 kappa_ip_state1 (HEAD rc 1,
  fix rc 0; reactions = Pro within 1e-4), test-2014 validation_3_kappa
  (HEAD rc 1, fix rc 0; Pro validation_3 displacements 1.25 % -> 0.06 %,
  kap node 1 11.6 % -> 2.7 %). build_safe TN_CHECK_ONLY=1 224 OK /
  0 FALLO, runs 226/226; audit_db_buffers exit 0. Corpus A/B wu1 vs
  wu2 (omen f4t6f_runs): rc changes = the 2 new tests 1 -> 0 and sfnet
  examp3 0 -> 1: post_node_result vely -126.431 (Pro on the translated
  input -126.374, 0.045 %; before -125.477) vs target -124.8 +- 1, a
  2014 GNU calibration -> owner decision (recalibrate as a stale sfnet
  target or not). Pro ground11_nonsaturated OK -> TIMEOUT was load
  (alone 96.9 s wu1 / 97.1 s wu2). BACKLOG C68 drafted (kapsh, cap1pc,
  rho, damage, Maxwell stresses: same mask omission).
  T6f WU2 resumed after the owner decision 2026-10-12 (recalibrate sfnet
  examp3 to the Professional value). `c3e6d2a` fix(plasti): kappa as an
  integration-point state (elem.cc masks + kap_indx, materi.cc IP write
  and nodal numerator with strain_recovery_weight outside the h guard,
  general.cc denominator); tests kappa_ip_state1 + validation_3_kappa
  (RED HEAD rc 1, GREEN rc 0). sfnet examp3 (local + omen copies,
  untracked): target -124.8 -> -126.374 +-1 with dated header (B16,
  CORPUS-STALE-TARGETS 8); GNU rc 0. A/B wu1 vs wu2 as recorded above
  (rc changes: 2 new tests 1 -> 0, examp3 0 -> 1 before recalibration).
  Clean build of c3e6d2a (build/tochnog, byte-identical to the A/B wu2
  binary): build_safe 224 OK / 0 FALLO, runs 226/226; audit_db_buffers
  exit 0. BACKLOG C68 open (kapsh, cap1pc, rho, damage, Maxwell msig).
  Next: T6f remaining items (mmc_tension, K14b, C7, K21, C57).
  T6f mmc_tension + validation_3 eptzz (2026-10-12, owner decision:
  class A, keep GNU, document), route: delegated writer (trigger: tests
  + 6 docs; no source change). `600a0b4` test(plasti): test-2026
  mmc_tension_vertex1/2 (test-2014 mmc_tension at dt 0.05 / 0.047;
  sigxx 46.19, sigzz 0 +- 0.05, eppzz 0 +- 1e-4; GNU sigzz -3.3e-10 /
  1.1e-9) and membrane_eptzz_total1/2 (vm1_ps at dt 0.01 / 0.1; eptzz
  analytic -0.0088249 +- 3e-4 / 2e-3, sigxx 27.4199 +- 0.05; GNU
  -0.0086300 / -0.0072102). Guard tests: rc 0 on HEAD by design, RED not
  applicable (GNU already correct); Professional values (sigzz 0.974 /
  0.367, eptzz -9.6e-5 / -7.9e-4) fail them. Finding: the analytic eptzz
  is -0.0088249, not -0.00863 (that is GNU at dt 0.01): GNU eppyy/eppzz
  split is first order in dt (dt 0.001 -0.0088050) -> BACKLOG C69.
  DIVERGENCES MC-VERTEX-ZZ + MEMBRANE-EPTZZ (A 20 -> 22, total 78);
  SEGUIMIENTO Descarte + register; new pages group_materi_plasti_mohr_coul
  and group_materi_membrane (user + developer, READMEs); BACKLOG C69,
  summary C64-C69, 47 open. build_safe TN_CHECK_ONLY=1 228 OK / 0 FALLO,
  runs 226/226. Next: T6f remaining items (K14b, C7, K21, C57).
  T6f C68 integration-point states (2026-10-12), route: delegated writer
  (trigger: 3 source files + 4 tests + docs; patch from a read-only
  investigation, scratch c68). The manual (4.10, 4.15, 4.17, 4.25, 4.27)
  treats kapsh, cap1 pc, damage, rho and the Maxwell stresses as point
  history like kappa -> GNU defect fix. `70235c6` fix(materi): elem.cc
  ELEMENT_DOF masks (feature-guarded), materi.cc one recovery-weight
  block after kappa outside the h guard + IP write (old h blocks
  removed), general.cc denominator; manual pages (damage, kapsh, cap1,
  Maxwell, kappa, nodal_strain_recovery; new materi_plasti_rho user +
  developer, READMEs). RED/GREEN on four-tria3 Pro probes: test-2026
  damage_ip_state1, kappa_shear_ip_state1, rho_ip_state1,
  maxwell_ip_state1 all rc 1 on c8737c8, rc 0 on the fix (vs Pro:
  damage reactions 1.9e-1 -> 5e-13, rho 1.1e-1 -> 7e-6, kapsh 1.3e-2 ->
  1.5e-4, Maxwell nodal msigyy 8.4e-2 -> 5e-4). cap1 manual-only (no Pro
  evidence). Corpus A/B omen f4t6f_c68_runs: no rc/class change, no
  node_dof change > 1e-10. Clean build; build_safe TN_CHECK_ONLY=1 232
  OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Docs commit:
  BACKLOG C68 resolved + C68a (Pro Maxwell vs analytic, Pro-suspect),
  C68b (cap1 model), C68c (MC h-s kapsh 1.2x), C68d (latent Maxwell
  index in the old-state read), 50 open; DIVERGENCES B 39 -> 42, total
  81; SEGUIMIENTO register. Next: T6f remaining items (K14b, C7, K21,
  C57).
  T6f C68d Maxwell old-state index (2026-10-12), route: delegated writer
  (trigger: source + test + build_safe + docs). Code reading confirmed
  only the old-state read loop of materi() had the pattern (the Maxwell
  loop overwrote indx before the rho/epi/epp reads). RED: test-2026
  maxwell_two_chains_state1 (two identical chains vs one chain with the
  summed stiffness; Pro 02-08-2026 equal to ~1e-12 on omen) rc 1 on
  30bb922 (node 1 reaction -13.32 vs -19.17); ASan/UBSan build of
  30bb922: no report (overrun inside dbl_array). `2327083`
  fix(materi): old_msig[m*MDIM*MDIM + indx], indx untouched; GREEN rc 0
  (1e-11). Clean build; build_safe TN_CHECK_ONLY=1 233 OK / 0 FALLO,
  runs 226/226; audit_db_buffers exit 0. Docs: BACKLOG C68d resolved,
  49 open; SEGUIMIENTO register + materi_maxwell_stress line; Maxwell
  developer page. Next: T6f remaining items (K14b, C7, K21, C57).
  T6f K14b + C7 (2026-10-12, owner decision: class A for both, keep GNU,
  document), route: delegated writer (trigger: 4 tests + build_safe + 10
  docs; no source change). Investigation: Engram 2742. `ae0691b`
  test(materi): test-2026 hypo_wolff_igs_oedometer1/2 (test-2014 hypo3
  at dt 1e-3 / 1e-5; sigxx converged reference -0.10466 from a Python
  integration of von Wolffersdorff + IGS, +- 2e-3 / 5e-4; GNU -0.103436 /
  -0.104235, dt 1e-5 runs in ~1 s; Pro -0.1629 / -0.10499),
  damage_mazars_eq_two_positive1 (hex8 diag(2e-4,-3e-4,1e-4), dam
  0.28610207515 +- 1e-8; GNU 0.286102; Pro 0.36409) and
  damage_mazars_vonmises1 (damage_stress1 bar + vonmises 3e6, dam
  0.925372225695 +- 1e-6, eppxx 1.87863e-4 +- 1e-8; Pro 0.96946 /
  1.368e-4). Guard tests: rc 0 on HEAD by design, RED not applicable (no
  behaviour change); sensitivity checked (a shifted target -> rc 1 with
  the actual value in tn.log). Finding: the Pro Mazars eps_eq s = 2 s +
  eps^2 is the 2014 `epseq += epseq + ...` defect the GNU fixed on
  2026-08-10 (manual-developer group_materi_damage_mazars "Bug report").
  DIVERGENCES K14b + C7 class A (A 22 -> 24, B 42 -> 41, total 82; K14 row
  keeps the hypo3 hyhis residual); SEGUIMIENTO Descarte x2 + register;
  BACKLOG C7 closed with the shared non-objective gamma note
  (`damage.cc:84-95`, no fix), 48 open (3 MEDIUM, 45 LOW); Differences
  sections of the hypo wolfersdorff / intergranularstrain and Mazars /
  damage effective stress pages (user + developer). build_safe
  TN_CHECK_ONLY=1 237 OK / 0 FALLO, runs 226/226; audit_db_buffers exit
  0. Next: T6f remaining items (K21, C57).
  T6f K21 WU1 nodal recovery (2026-10-12, owner-approved), route:
  delegated writer (trigger: 5 source files + 4 tests + build_safe + 10
  docs). Investigation Engram 2751. Probes (Pro 02-08-2026, operators
  fitted on element_intpnt_dof -> element_dof): distorted hex8 2x2x2 and
  curved quad9 3x3 are ALSO physical-coordinate fits (not implemented,
  BACKLOG C70); elastic axisymmetric quad4 confirms the r=0 rule (1.5e-13).
  RED on 7eed0db: test-2026 recovery_examp1_quad4_point_mean1 (node 19
  sigxx -0.256815 vs -0.2350738), recovery_incnav4_quad9_minimal1 (sigxy
  0.02/0.08 vs 0.05), recovery_tube1_axis1 (node 3 -66.40 vs -159.5805),
  recovery_quad4_distorted_plane1 (node 21 1.1829 vs 2.0014257), all rc 1.
  First build (ungated) failed test-2014 cd_method (convected
  materi_strain_total + derivatives, sigxy 0.744 -> 3.64, target 1+-0.3);
  a stress_recovery_absolute gate fixed it but also excluded axisym2
  (materi_strain_total); final gate stress_recovery_lagrange (convected
  stress only): cd_method unchanged, axisym_quad4_ring1 sigzz 10759.2155 =
  Pro. A per-call coordinate read cost +12 % on ground11_nonsaturated;
  the cache now re-reads only on the (inol 0, ipoint 0) call (+3 %).
  `2305659` fix(recovery): GREEN rc 0, whole stress field = Pro within
  1e-7. Clean build; build_safe TN_CHECK_ONLY=1 241 OK / 0 FALLO, runs
  226/226; audit_db_buffers exit 0. Corpus A/B omen f4t6g_rec_runs (head
  vs fin): rc changes only the 4 new tests (alternate1 timeout flake);
  closer 9, farther 0; no-Pro changes taylor3, disk1, examp18, axisym4,
  examp3, examp14, cd_geom (targets met); wall +0.9 %. Docs: BACKLOG C42,
  C64 resolved, C70, C71 new (open count unchanged); DIVERGENCES C42
  removed (closed by following the Pro), K21 re-scoped, C70/C71 added (B
  41 -> 42, total 83); SEGUIMIENTO register; AUDIT K21; user/developer
  nodal_stress_recovery, nodal_strain_recovery,
  integration_quadrature_and_stress_recovery. Next: K21 validation_1
  (stokes default, convection), C70, C57.
  T6f K21 WU1b gate removal + test rewrite (2026-10-12, owner-approved),
  route: delegated writer (trigger: 4 source files + 3 inputs + 8 docs).
  Investigation Engram 2756: `-ra 1 4 -ra` is a node LIST, cd_method /
  cd_geom rotated rigidly (physical stress 0, Pro cd_method 0 from t 0.6),
  old targets met by chance. `e3af9f2` test(materi): vely on every node +
  tension_direct 1.e10; GNU HEAD = new binary: cd_method 0.9999956 (1e-5,
  phi target kept), cd_geom 0 / 1 (1e-6), materi_direct_mc 1 (1e-6), rc 0
  (no RED for the rewrite: behaviour already correct, as expected). Pro
  02-08-2026 -updated_linear: 1 / 0,1 / 1 exactly; Pro -updated 18.04 /
  18.69 = GNU -updated, so the investigation's "Pro 18.0" was a memory
  translation artifact, not an MC-direct divergence (C17 note). RED for
  the gate removal: original cd_method rc 0 on HEAD (0.744), rc 1
  ungated (3.642). `e2b044f` refactor(recovery): gate removed (miscel.cc,
  general.cc, tochnog.h/-mod.h; materi.cc had no reference). Clean build;
  build_safe TN_CHECK_ONLY=1 241 OK / 0 FALLO, runs 226/226;
  audit_db_buffers exit 0. Corpus A/B omen f4t6g_ungate_runs (fin vs
  ungate): rc change only old cd_method 0 -> 1; cmat_gate (same ill-posed
  convected model) sig 5.9e-2 rel, target met, vs Pro 0.107 -> 0.118;
  Pro other / sfnet / test-2026 dbs unchanged (alternate1,
  cylinder_buckle near-120 s timeout flakes). Docs: BACKLOG C72 new
  (convected stress path, class B), C42/C17 notes; DIVERGENCES C72 (B 42
  -> 43, total 84); SEGUIMIENTO register; user/developer
  nodal_stress_recovery, integration_quadrature_and_stress_recovery.
  Next: C72 (convect IP stress, rewrite cmat_gate), K21 validation_1,
  C70, C57.
  T6f ranges (2026-10-12, owner-approved; audit Engram 2758): route
  delegated direct (writer). `35025dc` fix(input): range.cc range_scan
  (-STEP sign, from==to, several -from/-to blocks, reads through
  input_read_string, empty range error), range_expand value-block
  validation (descending -> Pro message, step <= 0, missing -to),
  input.cc non-integer inside a value range = error. RED on HEAD: 4 index
  tests rc 1, 2 value error tests rc 0 (silent); GREEN all 8; Pro
  02-08-2026 gives the same targets / errors (Pro rejects arithmetic
  words inside an index range, defines work). Audit probes rerun: all
  forms = Pro. Clean build; build_safe TN_CHECK_ONLY=1 249 OK / 0 FALLO,
  runs 226/226; audit_db_buffers exit 0. Corpus A/B omen f4t6g_range_runs
  (ungate vs range): no rc/class change, .dbs identical except
  bicg_csr_threads1 thread round-off. Docs: user/developer ranges (+
  READMEs), SEGUIMIENTO register, DIVERGENCES INT-DOT (class C 16 -> 17,
  total 85: 9.0 for an integer value outside a range, Pro rejects),
  BACKLOG C56 range.cc residual closed. Stale: manual-developer
  input_word_length.md line 39 still says range.cc pending (outside the
  edit surface). Next: C72, K21 validation_1, C70, C57.
  T6f K21 WU2 group_materi_stokes (2026-10-12, owner-approved), route:
  delegated writer (trigger: general.cc + 3 tests + cmat_gate + wiring +
  10 docs). Investigation Engram 2751 / 2765. Pro 02-08-2026 probes:
  absent record = -yes (validation_1 flow velx 0.1465713); the Pro reads
  index 0 for every group (-no/-yes, -yes/-no, record only on group 1,
  element 0 moved to group 1: always the index-0 value, 1e-12); Stokes vs
  convection -no residual = convected IP stress (density 0/1/2 no effect;
  GNU nodal stress convected but no velocity feedback) -> C72. First
  build (default only, index 0) hit the STOP RULE: test-2014 cmat_gate
  sigxy 1.532 -> 0.921 vs target 3.20774 +- 2 (Pro 0; ill-posed
  `-ra 1 4 -ra` rigid rotation, met by chance); corpus A/B omen
  f4t6g_stokes_runs (range vs stokes): rc change only cmat_gate, cd_method
  ept 4.6e-2 (target met), Pro other / sfnet / test-2026 .dbs unchanged.
  Owner decisions 2026-10-12: per element group read (manual; Pro index 0
  = class A STOKES-IDX), cmat_gate rewritten. `58bba14` fix(convection):
  general.cc default -YES + db(GROUP_MATERI_STOKES, element_group); tests
  RED on HEAD rc 1 / GREEN rc 0: stokes_default1 (Pro 0.1465713 +- 0.003,
  GNU 0.144159, HEAD -0.0570386), stokes_group1 (0.124105, Pro 0.1465713
  fails), stokes_group2 (all in a -no group = one-group -no -0.0570386,
  HEAD 0.144159, Pro 0.1465713 fails); cmat_gate simple shear with gate
  sigxy 250.5 / sigxx -49.5 (analytic = Pro exactly; no gate 150.5 /
  -149.5), GNU 250.4989 / -49.4998 tol 0.01, HEAD and new identical. No
  corpus input has an index != 0 stokes record, so the per-group build
  cannot change the A/B (not re-run). Clean build; build_safe
  TN_CHECK_ONLY=1 252 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0.
  Docs: user/developer group_materi_stokes (+ READMEs), SEGUIMIENTO line +
  register + Descarte, DIVERGENCES STOKES-IDX (A 24 -> 25, total 86), K21
  re-scoped, C72 widened; BACKLOG C72 (no new item); AUDIT K21. Next: K21
  WU3 validation_1 convection transient, C72, C70, C57.
  T6f C72 (2026-10-13, owner-approved; investigation Engram 2775), route:
  delegated writer. Premise corrected: neither the Professional 02-08-2026
  nor the GNU convects the integration-point stress (pulse bit-frozen).
  GNU bug fixed `c1a1ae9`: the convected nodal stress dof was incremental
  with a lumped mass weight F/volume (sri_stress_recovery_weight) against
  a convective weight h, so it moved at about v*volume_ip; now
  `stress_recovery_absolute = !materi_strain_total` and general() skips
  the convective terms of an absolute stress dof (option chosen: always
  absolute recovery, = Pro node stress; the weight-consistent transport
  was rejected because a transported nodal stress would still differ
  from the Pro and from the frozen IP state). Tests test-2026:
  `convected_stress_recovery1` RED on HEAD rc 1 (0.348732 / 0.904484 /
  0.6481) -> GREEN rc 0 (0.5 / 1.0 / 0.5 = Pro); `validation1_steady1`
  (t = 3 -0.0571954 vs t = 1 -0.0570386, tol 9.9e-4; HEAD identical, Pro
  +0.0034 fails). A/B on all `convection_apply -yes` inputs: velocities
  bit-identical, nodal stresses <= 1.3 % of max (validation_1, stokes_*);
  strconv1 byte-identical; cd_method / cmat_gate / materi_direct_wall
  unchanged. Clean build; build_safe TN_CHECK_ONLY=1 252 OK / 0 FALLO,
  runs 226/226 (new tests run manually rc 0, not yet registered:
  scripts/build_safe.sh outside the edit surface); audit_db_buffers exit
  0. Docs `cc245f7`: DIVERGENCES C72 rewritten + C72-V1 (class B 44,
  total 87), BACKLOG C72 nodal part resolved + C73 + C74, SEGUIMIENTO
  register, manual pages. Next: register the 2 tests in build_safe.sh,
  then C70, C57, C74.
  T6f C70 (2026-10-13), route: delegated writer (trigger: miscel.cc +
  7 tests + build_safe + 9 docs). Pro 02-08-2026 probes (distorted
  bricks, curved quad9/hex27/tet10, element_intpnt_coord ->
  element_dof): quad9 3x3, hex8 2x2x2, hex27 3x3x3 / 2x2x2, tet10
  5-point, prism6 6-point = linear least squares 1,x,y(,z) in PHYSICAL
  coordinates (residual <= 1e-12 of the scale); hex8 1 point = point
  value; axisymmetric quad9 3x3 / 2x2 with a node at r=0 = sum h sig /
  sum h (GNU 2x2 used the plain fit there); node = arithmetic mean.
  `39d0138` fix(recovery): recovery_fit_kind/shape/build/factor
  generalised to 3D, the F4 T9 natural table removed. Tests test-2026
  recovery_hex8_distorted1 (node 11 sigxx HEAD 1.1216 -> 1.3071691 =
  Pro), recovery_hex8_point_mean1 (2.9277 -> 2.3571531),
  recovery_hex27_distorted1 (3.4134 -> 3.6670886),
  recovery_tet10_curved1 (2.5316 -> 2.6187305), recovery_prism6_linear1
  (1.0575 -> 0.55 exact), recovery_quad9_curved1 (2.2832 -> 2.4332610),
  recovery_quad9_axis1 (-0.1123 -> 0.9518217): RED rc 1 on 421a3f4,
  GREEN rc 0, whole field = Pro <= 2e-12; registered in build_safe.sh.
  Clean build; build_safe TN_CHECK_ONLY=1 261 OK / 0 FALLO, runs
  226/226; audit_db_buffers exit 0. Corpus A/B omen f4t6g_c70_runs
  (421a3f4 vs 39d0138): no rc/classification change; force15 0.25 ->
  5e-11 and force9 3.4e-2 -> 1e-7 vs Pro, force11 hex27 nodal stress =
  Pro 7e-7 (C57 recovery part resolved; post_point part open), none
  farther; large1 +1.5 % wall. Found C75: the Pro prism6 points are the
  triangle corners x 2 Gauss levels, GNU interior (2/3,1/6,1/6) (nodal
  sigzz 1.6 % of the scale on a distorted prism). Docs `170e573`:
  DIVERGENCES C70 closed, C75 added (class B 44, total 87), C57
  re-scoped; BACKLOG C70 resolved, C75 new, C57/C67 updates; SEGUIMIENTO
  register; user/developer nodal_stress_recovery,
  integration_quadrature_and_stress_recovery, control_mesh_extrude.
  Next: C57 post_point, C75, C74, C71, K21 validation_1.
  T6f C75 (2026-10-13), route: delegated writer (trigger: Pro/GNU
  probes + 7 docs). Pro 02-08-2026 prism6 rule measured: triangle
  corners (1,0),(0,1),(0,0) x zeta 1/2 +- 1/(2 sqrt 3), weight 1/12, for
  the default, -normal, -maximal and inertia_apply -yes; -minimal = 2
  centroid points (GNU ignores the record for prism6 -> C76). Weights
  confirmed: a 20x2x2 prism6 cantilever's Pro nodal displacements equal a
  corner-rule FE model to every printed digit (GNU HEAD = interior rule).
  Manual 6.622 silent on the triangle points of a prism. Gate: both rules
  full rank and patch-test exact, but the corner rule doubles the bending
  error (cantilever vs Timoshenko, 20x2x2 / 40x2x4: -26.0 / -9.95 %
  against GNU -13.1 / -4.85 %; tip reaction +35 % vs +15 %), so the
  writer stopped before code. Owner decision 2026-10-13: keep the GNU
  rule, C75 = deliberate divergence class A. Docs `4ab9549`: DIVERGENCES
  C75 -> class A, C76 added class B (A 26, B 44, total 88); BACKLOG C75
  decided, C76 new; SEGUIMIENTO Descarte + register rows; user/developer
  integration_quadrature_and_stress_recovery and nodal_stress_recovery.
  No code, test or build change (suite not re-run). Next: C57 post_point,
  C74, C71, K21 validation_1, C76.
