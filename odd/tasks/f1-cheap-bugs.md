# F1 — cheap, isolated GNU-suspect bugs

## Objective
Close the cheap GNU-suspect items of the all-variables audit and the print
loop follow-up, each verified against the Professional 02-08-2026 or an
analytic value.

## Tasks
- [x] T1 C24: `print_fr.cc` (~229, print_frd) and `print_gm.cc` (~442, print_gmsh) loop `ipuknwn<nuknwn_` like the old print_dof (C23, `f22ba98`): repeated blocks + out-of-record reads. Fix like print_vt.cc / f22ba98; compare FRD/gmsh layout with Pro.
- [ ] T2 K03: default `convection_apply` (GNU `-yes`, general.cc:54; Pro manual 6.395 `-no`; `condif2` GNU T=1 vs Pro 0.5).
- [ ] T3 K04: Pro analytic targets the GNU misses (AUDIT §K04 table): condif7, force18, mohr_coul_direct8, mohr_coul_direct4, stress_pressure_factor, incnav3, tet10_2, ground9, axisym9, distri2 (elasti7 and tension_direct1 already fixed). Investigate one by one; fix the isolated ones, escalate the rest.
- [ ] T4 K13: direct/visco plasticity state output (`visc_pl4`: GNU eppxx 0 and f −1e10, Pro eppxx = ε − σ/E).
- [ ] T5 K02/K11: one-iteration lag of non-principal dofs (`spring4`) and nodal recovery of total strain (`tslv_part`).

## Rules
Principle: manual/analytic decides; Pro is the oracle when the manual is
silent; divergences need a user decision. Each task: RED/GREEN test,
`TN_CHECK_ONLY=1 ./scripts/build_safe.sh` green (226 runs, 101+ checks),
corpus A/B (test-2014 233/1, test-2026, sfnet 187/10), docs + SEGUIMIENTO.

## Progress
- Branch `dev/f1-cheap-bugs` from `055fbc7`.
- T1 (route: delegated direct writer, 2 non-trivial source files + test + docs): `print_frd()`/`print_gmsh()` loop over primary unknowns; FRD 2D vector 3rd component 0; gmsh names from dof type. Commit `ff9aa3e`. Evidence: new test `test-2026/print_frd_gmsh_dof_loop1` RED (old: 68 FRD blocks / 324 gmsh node_ fields, build_safe check FALLO) -> GREEN (5 / 17, OK); values = first block of each dof of Pro 02-08-2026; build_safe check-only 226/226 runs, 103 OK / 0 FALLO, rc 0; corpus A/B old vs new rc + .dbs identical (vs 233/10 = print_mesh_dof + 9 expected-error, t26 46/2 expected-error, sfnet 187/10). Remaining layout differences -> BACKLOG C27.
- T2 K03 (route: delegated direct writer): BLOCKED on a user decision, code NOT committed (left in the working tree). Verified: Pro manual 6.395 says default `-no`; Pro binary agrees (condif2 converted: no record 0.5, `-no` 0.5, `-yes` 1). But in the Pro `convection_apply` only governs MATERIAL convection wrt the mesh: `group_condif_flow` (condif1) and groundflow-velocity (condif5) heat convection stay on with `-no`/`-yes`/absent. Draft: `convection_apply_switch()` in general.cc (default -no, single source for general/beam/truss/spring/conspr) + condif_flow/gvel convection no longer gated by convection_apply. A/B vs T1: sfnet 187/10 -> 182/15 (condif2, condif3, examp1, incnav5, strconv1 lose their targets; condif2/3, examp1, incnav5 closer to Pro, strconv1 Pro RUNFAIL) + 12 output-only changes; validation-suite 233 -> 230 (cd_method, cmat_gate, materi_direct_wall: targets calibrated on the GNU default; cmat_gate closer to Pro, the other two Pro PARSE) + ~57 output-only changes (most round-off, some 1-13 %); test-2026 new convection_apply_default1/yes1 pass. Needs user decision: accept the Pro default and retarget / add explicit `-yes` to the GNU-corpus inputs, or keep the GNU default as a documented divergence.
