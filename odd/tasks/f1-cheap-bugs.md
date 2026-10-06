# F1 — cheap, isolated GNU-suspect bugs

## Objective
Close the cheap GNU-suspect items of the all-variables audit and the print
loop follow-up, each verified against the Professional 02-08-2026 or an
analytic value.

## Tasks
- [x] T1 C24: `print_fr.cc` (~229, print_frd) and `print_gm.cc` (~442, print_gmsh) loop `ipuknwn<nuknwn_` like the old print_dof (C23, `f22ba98`): repeated blocks + out-of-record reads. Fix like print_vt.cc / f22ba98; compare FRD/gmsh layout with Pro.
- [x] T2 K03: default `convection_apply` (GNU `-yes`, general.cc:54; Pro manual 6.395 `-no`; `condif2` GNU T=1 vs Pro 0.5).
- [ ] T3 K04 (5 fixed, 3 not GNU defects / already fixed, 2 parked for decision: mohr_coul_direct4 C30, tet10_2 C31): Pro analytic targets the GNU misses (AUDIT §K04 table): condif7, force18, mohr_coul_direct8, mohr_coul_direct4, stress_pressure_factor, incnav3, tet10_2, ground9, axisym9, distri2 (elasti7 and tension_direct1 already fixed). Investigate one by one; fix the isolated ones, escalate the rest.
- [ ] T4 K13: direct/visco plasticity state output (`visc_pl4`: GNU eppxx 0 and f −1e10, Pro eppxx = ε − σ/E).
- [ ] T6 C28: `node_force` sign/doubling (1D probe Pro +1, GNU −2), sibling of force18 (`fb2cf78`).
- [x] T7 C30 docs: record `mohr_coul_direct4` as a deliberate divergence (GNU follows manual 6.632).
- [ ] T5 K02/K11: one-iteration lag of non-principal dofs (`spring4`) and nodal recovery of total strain (`tslv_part`).

## Decisions
- 2026-10-06 user: C30 deliberate divergence (GNU = manual 6.632); C31 (tet10 nodal quadrature) and K20 (axisymmetric 1-point rule) moved to F4 (design + corpus A/B).

## Rules
Principle: manual/analytic decides; Pro is the oracle when the manual is
silent; divergences need a user decision. Each task: RED/GREEN test,
`TN_CHECK_ONLY=1 ./scripts/build_safe.sh` green (226 runs, 101+ checks),
corpus A/B (test-2014 233/1, test-2026, sfnet 187/10), docs + SEGUIMIENTO.

## Progress
- Branch `dev/f1-cheap-bugs` from `055fbc7`.
- T1 (route: delegated direct writer, 2 non-trivial source files + test + docs): `print_frd()`/`print_gmsh()` loop over primary unknowns; FRD 2D vector 3rd component 0; gmsh names from dof type. Commit `ff9aa3e`. Evidence: new test `test-2026/print_frd_gmsh_dof_loop1` RED (old: 68 FRD blocks / 324 gmsh node_ fields, build_safe check FALLO) -> GREEN (5 / 17, OK); values = first block of each dof of Pro 02-08-2026; build_safe check-only 226/226 runs, 103 OK / 0 FALLO, rc 0; corpus A/B old vs new rc + .dbs identical (vs 233/10 = print_mesh_dof + 9 expected-error, t26 46/2 expected-error, sfnet 187/10). Remaining layout differences -> BACKLOG C27.
- T2 K03 (route: delegated direct writer; user decision 2026-10-06: follow the Professional): default `convection_apply -no` via `convection_apply_switch()` (general.cc; readers general/beam/truss/spring/conspr), condif_flow/groundflow-velocity heat convection not gated by the switch (Pro scope, measured condif1/condif5). Commit `f745ce1`. Evidence: Pro manual 6.395 and Pro binary agree (condif2 0.5 absent/-no, 1 -yes); RED `convection_apply_default1` old temp 1 -> GREEN 0.5 = Pro; `convection_apply_yes1` 1 = Pro. test-2014 `cmat_gate`, `cd_method`, `materi_direct_wall` got explicit `convection_apply -yes` (targets calibrated on the old default; .dbs identical to T1 apart from the record). build_safe check-only 226/226, 105 OK / 0 FALLO, rc 0. A/B vs T1: vs 233/10 per-test unchanged (~57 output-only changes, Pro-closer where comparable); t26 48/2 all pass; sfnet 187/10 -> 182/15 (condif2/3, examp1, incnav5 = Pro; strconv1 Pro RUNFAIL) -> BACKLOG B10; 12 sfnet output-only changes.
- T3 K04 (route: inline parent + scratch builds; one work-unit commit per fixed item). Status re-run 2026-10-06 on build ba3e967 vs Pro 02-08-2026: `mohr_coul_direct8` already PASS (fixed by C14/C16); the other 9 still fail.
  - [x] condif7: `area.cc` `side_node_restriction_ok()` — the `_node`/`_element_node` restriction selects whole sides (all side nodes listed), every edge family. Evidence: condif7 20 -> 10 = Pro; probes heat/groundflow/convection/force_edge/force_edge_normal/support = Pro; tests `edge_node_restriction_heat1`/`_force1` RED->GREEN; build_safe check-only 226/226, 107 OK / 0 FALLO, rc 0; corpus A/B identical (vs 232/11, t26 50/2, sfnet 181/16 with this runner). Commit `4d36a89`.
  - [x] force18: `force.cc` `force_point_calculate()` adds the force (`+=`) on the ndim velocity equations (was `-=` over every principal unknown). Evidence: force18 velx −0.8 -> +0.8 = Pro; test `force_point_sign1` RED (−2) -> GREEN (1.999999902, tol 1e-6, point_el); build_safe check-only 226/226, 108 OK / 0 FALLO, rc 0; corpus A/B identical (+1 new t26 test). `node_force` sign/doubling -> BACKLOG C28. Commit `fb2cf78`.
  - [x] ground9: `interface.cc` `interface_element()` — interface groundflow permeability and condif conductivity per unit length (pair weight `w_ip*iface_measure`), interface groundflow capacity `0.5*C*w_i*L` per node and implicit (needed: the weighted permeability alone made `groundflow_interface` diverge). Evidence: ground9 1/3 -> 0.5 = Pro; widths 2/3 and thermal analog = Pro; `groundflow_interface` transient (t = 1/3/10, capacity 0/1/2, 2 x tall) = Pro; tests `interface_flux_measure_ground1`/`_condif1`/`interface_capacity_transient1` RED->GREEN; build_safe check-only 226/226, 111 OK / 0 FALLO, rc 0; corpus A/B: only `groundflow_interface` node 6 pres 3.79e-4 -> 3.96e-4 (= Pro). `condif11` parse -> BACKLOG C29. Commit `8fead02`.
  - [x] stress_pressure_factor: `stress.cc` `set_stress()` — trial dp factored when the step starts below the history, factor only strictly below `sph*(1-1e-10)` (rule measured on the Pro at 9 dt/factor combinations; one truncated-last-step Pro quirk documented). Evidence: −3.10 -> −3.05 = Pro; tests `stress_pressure_factor_reload1`/`_reload2` RED->GREEN; build_safe check-only 226/226, 113 OK / 0 FALLO, rc 0; corpus A/B identical. Commit `4580991`.
  - [x] incnav3: `general.cc` `general()` — MATERI density reduced by (1 − beta T) of `group_materi_expansion_volume` for the inertia (as materi.cc for gravity). Evidence: tvely −2.4976 -> −4.990415335 = Pro; test `expansion_volume_inertia1` RED->GREEN; build_safe check-only 226/226, 114 OK / 0 FALLO, rc 0; corpus A/B: sfnet `incnav3` legacy target −2.5 now fails (stale; Pro −4.9904 on it) -> BACKLOG B11, nothing else changes. Commit `7443ba8`.
  - [x] mohr_coul_direct8: already PASS on the re-run (C14/C16) — no action.
  - [x] mohr_coul_direct4: DECIDED 2026-10-06 deliberate divergence (T7). GNU = manual 6.632 (immediate return to c + |Fn| tan phi = 1.2027); Pro stays at 2.0 "plastic" (t 1..5), drifts to 1.095 (t 50) and 0.41 "elastic" (t 100). PRO-suspect -> BACKLOG C30.
  - [ ] tet10_2: PARKED, needs design. `polynom.cc` TET10 nodal rule is degree 1; scratch experiment with the degree-2 nodal rule (corners -1/20, mids 1/5) gives -1.0000000058 (exact), but negative corner weights break lumped inertia -> BACKLOG C31 with options.
  - [x] axisym9: not a GNU defect in K04 scope: uniform field, 0.16 % ring-force deficit from the forced 1-point axisymmetric quad4 rule (`polynom.cc`, audit K20) -> BACKLOG C32 (decision belongs to K20).
  - [x] distri2: not a defect: random-generator dependent target (both valid N(1,2e-3)) -> BACKLOG C32.
  - T3 close-out: final binary (clean build of `7443ba8` code) = build/tochnog; K04 re-run: condif7, force18, mohr_coul_direct8, stress_pressure_factor, incnav3, ground9 PASS; mohr_coul_direct4, tet10_2, axisym9, distri2 FAIL (parked / not defects, see above). build_safe check-only 226/226, 114 OK / 0 FALLO, rc 0. Corpus (audit runner, 60 s, vs base `ba3e967`): validation-suite 232/11 unchanged (only `groundflow_interface` node 6 pres -> Pro value), test-2026 48/2 -> 57/2 (+9 new tests), sfnet 181/16 -> 180/17 (`incnav3` stale legacy target, B11). Docs record `5420226`.
- T7 C30 (route: inline parent, docs only, no code): BACKLOG C30 -> deliberate divergence (decision + severity), C31/C32(K20) -> F4 (design + corpus A/B); AUDIT K04/K20 rows and recommendation updated; 'Differences with the Professional version' in `manual-user/` and `manual-developer/group_materi_plasti_mohr_coul_direct.md`; SEGUIMIENTO Descarte row. Evidence: structural readback of the diff (docs only, no build needed). Commit: see git log `docs(mohr_coul_direct)`.
