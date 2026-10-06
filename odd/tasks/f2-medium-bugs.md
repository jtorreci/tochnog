# F2 — medium GNU-suspect bugs (contact, truss/beam/support, power-law elasticity)

## Objective
Close the medium-size GNU-suspect clusters of the all-variables audit, each
verified against the Professional 02-08-2026 or an analytic value.

## Tasks
- [x] T1 K06: power-law elasticity semantics (`myoung6_e2`: GNU E=3000 → −4.0385 vs Pro E=1000+3000 → −5.3846; `myoung6_e3`; `mpower` poisson_power). Pro adds `group_materi_elasti_young` to the power law.
  - [x] T1a `poisson_power` replaces nu (stress.cc clears C before the rebuild) — `a270264`.
  - [x] T1b `young_power` + `young` additive in Pro — owner decision 2026-10-06: deliberate divergence (GNU = manual 2.2.2/6.661), no code change; docs `5f7c10c`.
- [ ] T2 K17: truss / beam / support_edge_normal divergence on Pro inputs (`trubea2/3` disz 2e5 vs 0.00333; `support_edge_normal_3` disy 124.75 vs 5.01e-4; `truss8`; `dynamic3`; `truss2/3/14`, `update2` ERRDUMP).
  - [x] T2a `support_edge_normal_3`: damping/density missing from the iteration matrix (area.cc) — `0b9ca9a`.
  - [x] T2b `truss14` / `truss2/3`: truss `-total` linear + no geometric stiffness (truss.cc) — `40567ff`.
  - [x] T2d `truss8` (C39): owner decision 2026-10-06 follow the Pro; nodes whose elements are all inactive by group_time are not driven by bounda_dof (bounda.cc) — `4ceed88`.
  - [ ] T2c PARKED: `trubea2/3` planar GNU beam vs Pro 3D beam (design, BACKLOG C37); `truss8` death-velocity reset (C39 residual), `update2` (not a truss, -updated robustness), `truss13` (truss_distribute_mpc not implemented), zero-force first iteration (C12) → C39; `dynamic3` → C11; support `_stiffness_freeze` → C38.
- [x] T3 K05 + C5: contact (`contac1` no interaction, GNU node 2 velx 0.03 vs Pro 0.02/0.01; `validation_10` sigyy ≈0 vs −21031; `contact_target_geometry` alias aborts "Data is allocated in a parallel loop" in 6 runs and drops an index-0 contact in examp9).
  - [x] T3a `contact_target_geometry`/`_switch` alias: db PUT inside the parallel contact loop (index 0 abort) + index > 0 dropped (contact.cc) — `3cef926`.
  - [x] T3b element-target contact (`contact_target_element_group`) without node_boundary: no mesh boundary generated → no interaction (`contac1/4`, `axisym3`) (area.cc + top.cc) — `50683a2`. Residuals → BACKLOG C40.

## Rules
Principle: manual/analytic decides; Pro is the oracle when the manual is
silent; divergences and many-target changes need a user decision. Each fix:
RED/GREEN test in test-2026 registered in build_safe.sh,
`TN_CHECK_ONLY=1 ./scripts/build_safe.sh` green (226 runs, 121+ OK, 0 FALLO),
corpus A/B with the same runner (validation-suite 233/10, test-2026, sfnet),
docs + SEGUIMIENTO + AUDIT/BACKLOG. Time-box ~1.5 h per root cause.

## Progress
- Branch `dev/f2-medium-bugs` from `e693018`.

### T1 (K06) — 2026-10-06, route: inline (one file, stress.cc, plus tests/docs)
- Re-run before (GNU build = HEAD d0fbe86, Pro 02-08-2026, converted inputs):
  `myoung6_e2` GNU −1.7308/−4.0385, Pro −2.3077/−5.3846; `myoung6_e3` GNU
  −0.4615/−1.0769, Pro −1.0385/−2.4231; `myoung6` GNU −0.7139/−1.6657, Pro
  −1.4351/−3.3485; `mpower` GNU −1.4389/−3.2924, Pro −0.4417/−1.4123.
- Without `group_materi_elasti_young`: Pro = GNU for all `myoung6*`; Pro
  history variable (calculated young) = 3000 on `myoung6_e2` while the
  stress needs 4000 → the Pro sums C(young)+C(E_law); manual 2.2.2/6.661
  say the law is the modulus → T1b decision (C36).
- `mpower`: the GNU summed C(young,0.3)+C(young,nu(p)) (forward-Euler replay
  reproduces −1.43887/−3.29236); Pro replaces nu (replay −0.43887/−1.41012).
  Fix stress.cc; GNU now −0.438872/−1.410125.
- Evidence: test-2026 `elasti_poisson_power_replace1` RED (old rc 1) /
  GREEN (new rc 0, Pro passes too); `mpower` retargeted −0.4389/−1.4101 ±0.02
  (only target change of T1). build_safe check-only (all F2 fixes): 226/226
  runs, 125 OK / 0 FALLO, rc 0. Corpus A/B (plain rc runner, old = HEAD
  build, new = all F2 fixes): validation-suite old 232/11 (mpower, new
  target) → new 233/10; test-2026 only the 4 new F2 tests change.

### T2 (K17) — 2026-10-06, route: inline per root cause (one source file each)
- Re-run before (GNU HEAD build vs Pro 02-08-2026, Pro inputs as distributed):
  `trubea2/3` GNU bicg breakdown, disz 2e5 (Pro 0.00333 exact);
  `support_edge_normal_3` GNU rc 1 disy 124.751 (Pro 5.00999e-4);
  `truss2/3` GNU solver breakdown, reaction velx −0.99997 (Pro 0.1/0.1);
  `truss14` GNU breakdown, disy 1e10 (Pro 0.70711/1.70711); `truss8` target
  passes both, only node 1 disx GNU 0.2 vs Pro 0.1; `dynamic3` target passes
  both, node 2 disx GNU 32.66 vs Pro 60.44; `truss13` GNU rc 1, node 9/10
  disx 0 (Pro 0.25/0.75); `update2` GNU diverges (Pro veliy 12.243).
- T2a support: bisection — without the damping step both agree (5.0e-4);
  without the freeze record Pro 5.00998e-4, GNU still 124.75 → damping.
  area.cc matrix carried only k*dt; c = 1e6 vs k*dt = 1e3 → overshoot.
  Fix: k*dt + c + d/dt. GNU now 5.00998004e-4 = Pro (no freeze); with the
  freeze record Pro 5.00999e-4 (freeze parse-only in GNU, BACKLOG C38).
  Test `support_edge_normal_damping_tangent1` RED rc 1 / GREEN rc 0 (Pro
  passes). Corpus: no rc change; tsup_* dumps identical.
- T2b truss: `truss14` (`-total`) fell into the linear branch (initial
  axis) → zero transverse stiffness; `truss2/3` (`-updated`) had no
  geometric stiffness. Fix: `-total` branch (current axis, EA/L0) + tension
  geometric stiffness N/L (I − n n) dt for -updated/-total. GNU truss14 =
  Pro (0.707107/1.707107, N 1.414214); truss2 reactions 0.1/0.1, veliy
  0.823578 (Pro 0.823577). Tests `truss_total_rotation1`,
  `truss_updated_transverse1` RED rc 1 / GREEN rc 0 (Pro passes both).
  Corpus: no rc change; dumps differ only by iteration path (bmom_truss
  3.6e-7; sfnet truss2/3 2.9e-5, truss7 3e-5). First iteration of a
  zero-force straight truss still prints "iterative solver broke down"
  (C12) and continues.

### Closure evidence (2026-10-06)
- Commits: `a270264` fix(stress) T1a, `0b9ca9a` fix(area) T2a, `40567ff`
  fix(truss) T2b; parking docs in the chore commit with this update.
- Final clean build of HEAD (`git archive`, `build_safe.sh --clean`) is
  byte-identical to the A/B "new" binary; copied to `build/tochnog`.
  `TN_CHECK_ONLY=1 ./scripts/build_safe.sh`: 226/226 runs, 125 OK / 0
  FALLO, rc 0 (baseline 121 OK + 4 new tests).
- Corpus A/B (same plain rc runner, timeout 60 s, old = d0fbe86 build):
  validation-suite 232/11 → 233/10 (mpower with the new target); test-2026
  64/6 → 68/2 (the 4 new tests; 2 expected-error remain); sfnet 181/16 →
  181/16 (no rc change); Pro `other` 231/99 → 234/96 (support_edge_normal_3,
  truss2, truss3 now pass; no regression).
- Open decision: T1b (C36). T3 not started (separate agent).

### User decisions applied (2026-10-06/07)
- T1b / C36: deliberate divergence, the GNU keeps law = modulus; code and
  `myoung6*` targets unchanged. BACKLOG C36 decided, SEGUIMIENTO
  "Descarte de features" row, "Differences with the Professional" in the
  young_power user/dev pages. Commit `5f7c10c` (docs).
- `mpower` retarget (-0.4389/-1.4101 +/- 0.02, `a270264`): APPROVED by
  the owner; recorded in the SEGUIMIENTO poisson_power row (`5f7c10c`).
- T2d / C39 truss8, route inline (bounda.cc + shared helper in group.cc,
  elem.cc reuses it): root cause `bounda()` applied `bounda_dof` to every
  node regardless of element activity, while `elem()` skipped the unborn
  truss -> node 1 driven 0.05/0.1 before birth. Fix: static
  `bounda_node_dormant()` (NULL when no group_time record) + skip of the
  prescribed-dof branch for dormant nodes. Pro truss8 node 1 disx
  0/0/0.05/0.1 = GNU after (before 0.05/0.1/0.15/0.2). Test
  `group_time_bounda_dormant1` RED rc 1 old / GREEN rc 0 new, Pro rc 0.
  Death probe (`group_time 0 0. 1.`): Pro zeroes the node velocities after
  death, GNU unchanged (keeps the old velocity) -> C39 residual, not done.
  Corpus A/B (plain rc runner, timeout 60 s, -P 6, old = b042128 build):
  validation-suite 233/10 -> 233/10; test-2026 68/3 -> 69/2 (only the
  new test); sfnet 181/16 -> 181/16; Pro `other` 237/93 -> 237/93 (truss8
  rc 0 both). Stdout identical except timeouts/peak-memory noise. Only
  truss8 and the new test contain `group_time` in all corpora. Commit
  `4ceed88`.
- Final clean build of `4ceed88` (`git archive`, `build_safe.sh
  --clean`) copied to `build/tochnog`; `TN_CHECK_ONLY=1
  ./scripts/build_safe.sh`: 226/226 runs, 126 OK / 0 FALLO, rc 0.
- Next: T3 (contact) not started.

### T3 (K05 + C5) — 2026-10-07, route: inline per root cause
- Re-run before (GNU = build/tochnog of `841f6f3`, Pro 02-08-2026): Pro
  dialect `contac1` GNU rc 1 node 2 velx 0.03 (Pro 0.02/0.01);
  `contac4`, `axisym3` GNU rc 1 (no contact); `contac2` GNU rc 1 x 1.02
  (Pro 1.0); `validation_10` GNU rc 1 sigyy ~0 (Pro -21031);
  `pressed_rubber_cylinder` max |dof - Pro| 0.10; renamed GNU inputs
  (`contact`, `contact_block`, `contact_heatgen`, `examp6/9/24`) abort
  "Program error ... contact_geometry. Data is allocated in a parallel
  loop". Two root causes: (a) the alias, (b) no node_boundary.
- T3a root cause: `contact.cc` `parallel_contact()` start (old lines
  71-86) copied `CONTACT_TARGET_GEOMETRY 0` with `db(...PUT)` while
  `parallel_active` (called via `parallel_sys_routine`, top.cc:560) →
  `db_allocate` abort; indices > 0 never read. Fix: read-only
  `contact_geometry_get()` (both names, every index, GNU name wins).
  After: validation_10 sigyy -21232 (Pro -21031, target passes);
  pressed_rubber_cylinder |dof - Pro| 6e-7; contac2 1.0 = Pro; renamed
  inputs byte-identical `node_dof` to the GNU-name originals.
  Tests `contact_target_geometry_index1/_index0/_switch1` RED rc 1 old
  (1.02 free flight / parallel abort) / GREEN rc 0 new; Pro passes all 3.
  build_safe check-only (RC1 binary): 226/226 runs, 129 OK / 0 FALLO,
  rc 0. Corpus A/B (plain rc runner `corpus.py`, timeout 60 s, 6
  workers, old = `841f6f3` build, new = T3a): validation-suite 233/10 →
  233/10; test-2026 69/2 → 69/2 (new tests not in the old run); sfnet
  181/16 → 181/16; Pro `other` 233/97 → 235/95 (`contac2`,
  `pressed_rubber_cylinder` rc 1 → 0; 18 timeouts identical both sides,
  this run was concurrent). dbs identical except those two (filtering
  time stamps and the pre-existing uninitialized `node_node` /
  `element_interface_intpnt_materi_tension_status` garbage).
- T3a commit `3cef926`.
- T3b root cause: `contact.cc` `parallel_contact()` loops element targets
  only when `db_max_index(NODE_BOUNDARY) >= 0` and uses a side only when
  its nodes are node_boundary; the GNU generates node_boundary only in
  macro/adjust/delete (or by hand in GNU 2014 inputs), the Professional
  computes it in every run. Fix: `area.cc` `contact_node_boundary_set()`
  (nodes of unshared sides, `element_normal_sides` table) called from
  `top.cc` `iteration_start()`, only with `contact_target_element_group`
  and no node_boundary record. node_boundary sets = Pro (contac1 1,2,3,5;
  2D 3x2 quad4 block all but the 2 interior nodes). contac1 velx
  0.02/0.01 = Pro (before 0.03/0); contac4 = Pro exactly; axisym3 velx
  -9.99981e-4 = Pro. Tests `contact_element_target_boundary1` (1D) and
  `_boundary2d1` (2D) RED rc 1 old / GREEN rc 0 new; Pro passes both.
  build_safe check-only (RC2 binary): 226/226 runs, 131 OK / 0 FALLO,
  rc 0. Corpus A/B (same runner, old = T3a build): validation-suite
  233/10 → 233/10; test-2026 69/2 → 74/2 (the 5 new T3 tests); sfnet
  181/16 → 181/16 (examp14 124 → 1 = load noise, stale target B2);
  Pro `other` 235/95 → 242/88: `axisym3`, `contac1`, `contac4` rc 1 → 0,
  plus `dynamic5`, `large1`, `large3`, `taylor3` 124 → 0 (timeouts of
  the concurrent old run; this run was alone). dbs identical except the
  three contact inputs. Net vs the F2 baseline (Pro other 237/93): +5
  (contac1/2/4, axisym3, pressed_rubber_cylinder), no regression.
- Not done (BACKLOG C40): Pro error for single-layer targets, hex8 sides
  (`number_of_side_boundary_nodes==ndim`), element_boundary record.
- T3b commit `50683a2`. Final clean build of `50683a2` (`git archive`,
  `build_safe.sh --clean`) is byte-identical to the A/B RC2 binary;
  copied to `build/tochnog`; check-only 226/226 runs, 131 OK / 0 FALLO.
- Next: F2 closed except the parked T2c items; no open T3 decision.

