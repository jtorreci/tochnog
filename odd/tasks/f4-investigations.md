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
- [x] T6 Open-cause clusters: K14 (hypoplasticity history/stress path), K15 (MC/DP/von Mises), K24 (3D consolidation `large2`), K21 (Navier-Stokes), C7 (Mazars), C11 (dynamics), C3 (3D section sign), C6 (iface_gen coverage). Investigations, one by one. State 2026-10-17: K14 slots done (K14a), `hypo3` stresses class A K14b, hyhis residuals open/deferred (AUDIT K14); K15 resolved or class A (`mohrcou2` parked); K24 class A; K21 recovery + stokes default resolved, `validation_1` decided C72-V1, convected-stress residuals = BACKLOG C72; C7 class A; C11 resolved; C3 sign resolved + square-face tie class A; C6 resolved. Remaining residuals tracked in BACKLOG/DIVERGENCES, no open T6 cluster.
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
- [x] T6f Rest by relevance (K15 `validation_14` done `1d446b3`; K15 validation_3 kap done `c3e6d2a`; C68 IP states kapsh/cap1/damage/rho/Maxwell done `70235c6`; C68d Maxwell old-state index done `2327083`; K21 recovery part done `2305659`, gate removed `e2b044f`; `-ra` ranges done `35025dc`): K15 `validation_14` (re-scoped by T6c: pass-1 hole-boundary nodal stress/plastic strain up to ~20 % at 8 nodes on the r=0.01 circle -> `refine_locally` picks a different element -> pass 2 on a different mesh), K15 mmc_tension and validation_3 eptzz closed 2026-10-12 as class A divergences (`600a0b4`), K14b hypo3 and C7/C7b Mazars closed 2026-10-12 as class A divergences (`ae0691b`), K21 (recovery done; validation_1 stokes default + convection remain), C57 (hex27 nodal stress resolved by C70; post_point open), C70 done `39d0138` (all Gauss rules + 3D physical fit), C71 done `bea7599` (on-axis Lobatto -> Gauss; Lobatto recovery weighting -> C79), C79 done `ab4ef78` (Lobatto recovery = Pro + msf on-axis replica; found C80 `-reduced`), C80 done `995de49` (`-reduced` = Pro; default factor 10 class A; next C80.1 selective/B-bar), C75 (prism6 points). Parked: K14c (ISA, class D), mohrcou2 sfnet (input artifact). State 2026-10-17: every listed item done or decided (K21 `validation_1` = C72-V1; C57 post_point done `cf37de3`; C75 class A); the LOW queue closed below; parked: K14c (ISA, class D), mohrcou2 sfnet (input artifact).
- [x] T6f C72 convected nodal stress + validation_1 steady guard (owner-approved 2026-10-13; route: delegated writer, trigger: general.cc + 3 tests + 9 docs). Done `c1a1ae9` (fix + tests), `cc245f7` (docs). Tests registered in `scripts/build_safe.sh` `421a3f4`.
- [x] T6f C70 nodal recovery of the other Gauss rules and the 3D elements in physical coordinates (route: delegated writer; trigger: multi-file code + 7 tests + build_safe + 9 docs). Done `39d0138` (fix + tests + registration), `170e573` (docs). Found C75 (prism6 integration points).
- [x] T6f C75 prism6 integration rule: investigation + owner decision + docs (route: delegated writer; trigger: Pro/GNU probes + 7 docs). Gate stop before code; owner 2026-10-13: keep the GNU interior rule (class A). Done `4ab9549` (docs). Found C76 (`-minimal` ignored for prism6).
- [x] T6f C74 Bi-CG breakdown message on `validation_1` (route: delegated writer; trigger: Pro/GNU probes + corpus A/B + so.cc/so_bicg.cc + test + 5 docs). Investigation stop, owner 2026-10-13: option A (one warning per LU-rescued breakdown), option B -> C77. Done `cf07b19` (fix + test + build_safe), `e5f50fc` (docs).
- [x] T6f C57 `post_point` on the quad9 -> hex27 extruded mesh (route: delegated writer; trigger: Pro/GNU probes + corpus A/B + mesh.cc + test + 9 docs). Investigation stop (mesh.cc outside the surface), owner authorized 2026-10-13; curved-element rule C78 decided deliberate divergence class A. Done `cf37de3` (fix + test + build_safe), `30f2621` + `1f5e82c` (docs).
- [x] T6f C79 Lobatto nodal recovery + C71 `msf_element_rule()` replica (route: delegated writer; trigger: Pro/GNU probes + corpus A/B on omen + miscel.cc/polynom.cc/calcul_force.cc/headers + 6 tests + 9 docs). Done `ab4ef78` (fix + tests + registration), `6094739` (docs). Found C80 (`-reduced` refused by the GNU).
- [x] T6f C77 iterative-solver relative stop + BiCGStab (relative stop done `d567672` + docs `5091fe7`; BiCGStab vs Bi-CG choice returned to the owner) (owner-approved 2026-10-15: relative stop `||r||/||b||` + relative breakdown + `eps*||A||*||x||` floor, no stagnation window / best iterate; respect `bicg_error_minimum`; guarded direct fallback SuperLU -> band LU with a memory budget; BiCGStab (or GMRES(m)) replaces Bi-CG for the non-symmetric path; CG gets the relative stop only if measured safe; ILU/AMG -> F6; MKL PARDISO discarded). Route: delegated writer (trigger: so_bicg.cc + so.cc + tests + build_safe + 6 docs, corpus A/B on omen, local large-mesh timing). Subtasks:
  - [x] C77a RED test(s) under test-2026 (scaled validation_1, HEAD absolute floor accepts a far-from-converged solve; target SuperLU) + build_safe registration.
  - [x] C77b relative stop + relative breakdown + eps floor + `bicg_error_minimum` honoured + guarded fallback (SuperLU, then band LU within a memory budget, clear error otherwise); C74 warning adapted.
  - [x] C77c (owner decision 2026-10-15: keep Bi-CG with the relative stop; BiCGStab revisited with ILU in F6, BACKLOG C81) BiCGStab for the non-symmetric path; measure vs relative-stop Bi-CG (86 corpus cases + large C13/C23/C36, E15/E25); CG relative stop measured.
  - [x] C77d suite targets that move: verdict per target vs SuperLU/Pro; omen corpus A/B.
  - [x] C77e docs: user/dev solver pages, DIVERGENCES, SEGUIMIENTO (register + Descarte MKL PARDISO), BACKLOG (C77 closed, F6 ILU/AMG item).
- [x] T6f C80 `group_integration_method -reduced` + `group_integration_method_reduced_factor` (owner decisions 2026-10-16: default factor 10 per manual 6.621, class A divergence with Pro's 1; hourglass warning in the user manual; selective/B-bar -> new C80.1, not implemented here). Route: delegated writer (trigger: headers + database.cc + polynom.cc + miscel.cc + calcul_force.cc + delete.cc + 5 tests + build_safe + 7 docs, corpus A/B on omen). Done `995de49` (feat + tests + build_safe), `6c8269f` (docs). Subtasks:
  - [x] C80a RED tests under test-2026 (quad4 f=2 and f=10, hex8 f=2, quad9 `-reduced -minimal`: Pro targets 4e-6; no-record default f=10: GNU-consistent target) + build_safe registration.
  - [x] C80b enums REDUCED + GROUP_INTEGRATION_METHOD_REDUCED_FACTOR (both headers), database.cc name + record (factor > 0), pol() rule table, recovery kind 1 at the reduced points, msf/mc replicas, C71 on-axis switch for `-reduced`.
  - [x] C80c verification: clean build, build_safe TN_CHECK_ONLY, audit_db_buffers, omen corpus A/B.
  - [x] C80d docs: user/dev integration pages, DIVERGENCES (C80 closed + default-factor class A), SEGUIMIENTO (register + Descarte), BACKLOG (C80 closed + C80.1 OPEN next).
- [x] T6f C80.1 selective volumetric reduced integration (owner decision 2026-10-16: extend the GNU-only `group_element_selective_reduced_integration` with `-volumetric` (B-bar, volumetric part at 1 centroid point, deviatoric full Gauss) and `-both` (shear + volumetric at 1 point); `-yes` unchanged; opt-in only). Route: delegated writer (trigger: miscel.cc + materi.cc + database.cc + check.cc + calcul_force.cc + headers + 4 tests + build_safe + 6 docs, omen corpus A/B). Done `b47c39c` (feat + tests + build_safe), `7504e47` (docs). Subtasks:
  - [x] C80.1a RED tests under test-2026 (Lame cylinder plane strain nu 0.49 and 0.4999 `-volumetric`, coarse plane-strain cantilever `-both`, hex8 one-layer Lame `-volumetric`; analytic targets, tolerance from the measured discretisation error); `-yes` covered by the existing qsri/qsri3d tests; build_safe registration.
  - [x] C80.1b enums VOLUMETRIC + BOTH (both headers), database.cc names, sri_active() mode mask (shear 1, volumetric 2; volumetric dropped for membrane groups), materi() volumetric/deviatoric split (kappa = m^T D m / 9, centroid divergence row, D-c feedback, B-bar stored stress), calcul_force shear-only replacement for the shear bit, check.cc value check.
  - [x] C80.1c measurement vs analytic (Gauss, `-yes`, `-volumetric`, `-both`, Pro `-reduced` f=2/10 from the C80 study), hourglass single quad4 + 4x4 point-load block, element eigenvalues.
  - [x] C80.1d verification: clean build, build_safe TN_CHECK_ONLY, audit_db_buffers, omen corpus A/B.
  - [x] C80.1e docs: user/dev SRI pages (GNU-only note, comparison vs `-reduced`), SEGUIMIENTO register, DIVERGENCES (only if GNU-only extensions are listed there), BACKLOG C80.1 closed.
- [x] T6f C80.2 `-lambda` (Hughes 1977 lambda/mu split) + removal of `-both` (owner decision 2026-10-16: only lambda*m*m^T at the centroid point, the 2*mu part full Gauss, linear elastic; `-both` removed and rejected; `-yes` / `-volumetric` unchanged; gate: no material regression vs `-volumetric`/Gauss on the C80.1 problems, no spurious mode, before committing). Route: delegated writer (trigger: miscel.cc + materi.cc + database.cc + input.cc + headers + tests + build_safe + docs, omen corpus A/B). Done `5993499` (feat + tests + build_safe), `7866a53` (docs); gate passed. Subtasks:
  - [x] C80.2a RED tests under test-2026 (`-lambda` coarse plane-strain cantilever 10x1 nu 0.3, Lame cylinder nu 0.4999; `-both` expected-error test) + build_safe registration (the `-both` test removed).
  - [x] C80.2b enum BOTH -> LAMBDA in place (both headers, numbering unchanged), database.cc name, sri_active() SRI_LAMBDA bit (dropped on membrane groups), materi() lambda coefficient in the C80.1 centroid machinery (isotropic tangent check), input.cc clear `-both` message.
  - [x] C80.2c gate measurement (C80.1 problems + pure bending + hex8 cantilever + hourglass + eigenvalues + mesh convergence).
  - [x] C80.2d verification: clean build, build_safe TN_CHECK_ONLY, audit_db_buffers, omen corpus A/B.
  - [x] C80.2e docs: user/dev SRI pages (comparison table, guidance), SEGUIMIENTO register, BACKLOG (C80.2 note on C80.1).
- [x] T6f C67 nodal recovery of the bar2/bar3/bar4 and quad16 rules (hex64: rejected by both codes). Pro measured (65 probes, every node): Gauss and -reduced and Lobatto -minimal = physical linear fit `1,x(,y)`; Lobatto -maximal = point value of the same tensor index (also bar4/quad16, whose interior points are not at the nodes); arithmetic mean over elements. Route: delegated writer (trigger: miscel.cc + polynom.cc + tests + build_safe + docs, omen corpus A/B). Done `d875874` (fix + tests + build_safe), `aee381e` (docs). Subtasks:
  - [x] C67a RED tests under test-2026 (bar2 2-point Gauss, bar3 Gauss / Lobatto -minimal / -reduced, bar4 Gauss / Lobatto -maximal, quad16 Gauss / Lobatto -maximal / -minimal; Pro targets 4e-6) + build_safe registration.
  - [x] C67b miscel.cc recovery_fit_kind/recovery_fit_shape for 1D bars and the 4-node-per-direction elements; polynom.cc -reduced Lobatto rules of every tensor element recovered by the fit.
  - [x] C67c verification: clean build, build_safe TN_CHECK_ONLY, audit_db_buffers, every probe node = Pro, omen corpus A/B.
  - [x] C67d docs: user/dev nodal_stress_recovery, integration_quadrature_and_stress_recovery, group_integration_method_reduced; DIVERGENCES, SEGUIMIENTO register, BACKLOG C67.
- [x] T6f C67 remainder + C76 nodal recovery of the remaining rules on the `h` fallback (tet4 `-maximal`, prism15, tria3/tria6 axisymmetric variants, any other element/rule combination of `pol()` not covered by `sri_stress_recovery_weight()`) and the prism6 `-minimal` 2-point rule. Route: delegated writer (trigger: Pro/GNU probes + miscel.cc/polynom.cc (+ msf/mc replicas if a rule changes) + tests + build_safe + docs, omen corpus A/B). C75 (GNU interior prism6 rule, class A) stays. Done `dc4a59a` (fix + tests + build_safe), `5bd6f90` (docs). Open owner questions: C82 (prism15 `-minimal`), C83 (tria3/tet10 `-maximal`), C66 (tria6 rule + lumping). Subtasks:
  - [x] C76a enumerate the fallback combinations from `pol()`/`recovery_fit_kind()`; Pro probes (`scratchpad/c76/`) of every combination: points, weights, per-element operator, averaging; C76 points/weights vs GNU and the analytic error if the rule differs.
  - [x] C76b RED tests under test-2026 (Pro targets 4e-6) + build_safe registration.
  - [x] C76c implementation (miscel.cc recovery, polynom.cc rule if needed, replicas consistent).
  - [x] C76d verification: clean build, build_safe TN_CHECK_ONLY, audit_db_buffers, omen corpus A/B.
  - [x] C76e docs: user/dev nodal_stress_recovery (+ integration page), DIVERGENCES C76, SEGUIMIENTO register, BACKLOG C67/C76.
- [x] T6f C82 + C83 (owner decisions 2026-10-16: C82 prism15 `-minimal` keep 21 points, class A; C83 part 1 tria3/tet10 `-maximal` keep 3/5 points, class A divergence from the MANUAL (GNU = Pro); C83 part 2 axisymmetric tria3 `-maximal` follow the Professional: 3 nodal points, 1 point in an element with a node on the axis). Route: delegated writer (trigger: Pro/GNU probes + polynom.cc + test + build_safe + 6 docs, omen corpus A/B). Done `72f78bd` (fix + test + build_safe), `8e9606d` (docs). Subtasks:
  - [x] C83a why the GNU forced 1 point (blame/comments); Pro default (no record) for the axisymmetric tria3; replicas (`msf_element_rule()`, `mc_element_rule()`) tensor-only check; tet10 positive-weight 10-point rule check.
  - [x] C83b RED test under test-2026 (axisymmetric tria3 `-maximal`, elements on and off the axis, Pro targets 4e-6) + build_safe registration.
  - [x] C83c polynom.cc `pol()` axisymmetric tria3 `-maximal` = 3 nodal points, on-axis node -> 1 point (recovery unchanged if the probes agree).
  - [x] C83d verification: clean build, build_safe TN_CHECK_ONLY, audit_db_buffers, probes = Pro, omen corpus A/B.
  - [x] C83e docs: DIVERGENCES (C82, C83 part 1 class A; C83 out of class B), SEGUIMIENTO (Descarte C82 / C83 part 1, register C83 part 2), user/dev nodal_stress_recovery, integration pages, prism15; BACKLOG C82/C83 closed.
- [x] T6f C66 tria6 rule + lumping (owner 2026-10-16: MEASURE how the Professional lumps the mass before any lumping decision). Route: delegated writer (trigger: Pro/GNU mass probes + lumping-site inventory + possible polynom.cc/miscel.cc/general.cc + tests + docs, omen corpus A/B). Phase 2 gate: implement only if the Professional lumping is sound (positive, total-preserving); otherwise stop before any source edit. Subtasks:
  - [x] C66a GNU lumping-site inventory (general.cc inertia block, recovery weights, calcul_force/elem/interface lumped terms).
  - [x] C66b Pro lumped masses (point force, other nodes fixed, v = F dt / m): tria3, quad4, quad9, tet10, tria6 default/-maximal/-minimal, plane and axisymmetric; condif capacity cross-check.
  - [x] C66c verdict vs row-sum / HRZ / nodal quadrature / equal split -> gate result.
  - [x] C66d (only if the gate passes) implementation + RED/GREEN tests + verification: done in phase 2 below.
  - [x] C66e docs: BACKLOG C66 (measured lumping), DIVERGENCES C66 row, BACKLOG new item for the `post_element_force` tria3/tet path (C84).
- [x] T6f C66 phase 2 (owner decision 2026-10-16: manual 6.623 rule, 7 points default/`-maximal`, 4 points `-minimal`, physical linear-fit recovery + arithmetic mean = Professional; HRZ lumping of the h-weighted tria6 dofs, plane corner A/19 and mid-side 16A/57, axisymmetric r-weighted diagonal scaled to 2 pi A r_c, deliberate class A divergence TRIA6-HRZ; accept axisymmetric tria6, match the Professional on-axis behaviour). Route: delegated writer (trigger: polynom.cc + miscel.cc + general.cc + materi.cc + headers + tests + build_safe + docs, Pro probes, omen corpus A/B). Done `175faf1` (fix + tests + build_safe), `4b65b07` (docs). Subtasks:
  - [x] C66f Pro probes: tria6 point order/weights, axisymmetric on-axis tria6 (rule switch?), r=0 hoop term.
  - [x] C66g RED tests under test-2026 (plane/curved recovery 7 and 4 points, axisymmetric recovery, HRZ masses plane/axisymmetric; Pro / hand targets) + build_safe registration.
  - [x] C66h implementation: pol() rule + axisymmetric acceptance, recovery_fit_kind/shape, HRZ weight (general.cc inertia block + materi.cc state numerators), headers; replicas tensor-only check.
  - [x] C66i analytic HRZ sanity check (row sum vs HRZ vs analytic).
  - [x] C66j verification: clean build, build_safe TN_CHECK_ONLY, audit_db_buffers, omen corpus A/B of every tria6 case.
  - [x] C66k docs: user/dev nodal_stress_recovery, integration pages, TET10-HRZ/mass-lumping page, DIVERGENCES TRIA6-HRZ, SEGUIMIENTO register + Descarte, BACKLOG C66 closed.
- [x] T6f C85 + C82 consistency check (owner decisions 2026-10-17: C85 apply the tria6 HRZ weight (`h_lump` pattern) to the `groundfl.cc` groundflow-velocity numerator and the `condif.cc` lumped absorption, and switch the `general.cc` `groundflow_velocity` denominator to HRZ; C82 check: tria6 `-minimal` 4-point stiffness eigenvalues, stop and report if negative). Route: delegated writer (trigger: groundfl.cc + condif.cc + general.cc + 2 tests + build_safe + docs, eigenvalue study, omen corpus A/B). Done `50fd4c7` (fix + tests + build_safe), `eb00f27` (docs). C82 check STOPPED: tria6 `-minimal` is indefinite on curved / axisymmetric elements at high nu (owner question, no rule change, no contrast documented). Subtasks:
  - [x] C85a RED tests under test-2026 (tria6 groundflow velocity from a linear pressure; tria6 condif absorption decay; analytic / hand targets) + build_safe registration; RED on c6fa9aa.
  - [x] C85b implementation: groundfl.cc / condif.cc `h_lump`, general.cc exception dropped; GREEN; other elements byte-identical.
  - [x] C85c C82 consistency: tria6 4-point (-27/48) element stiffness eigenvalues (straight/curved, nu 0.3/0.49, plane strain/axisymmetric): zero modes and negative eigenvalues, contrast with prism15.
  - [x] C85d verification: clean build, build_safe TN_CHECK_ONLY, audit_db_buffers, omen corpus A/B.
  - [x] C85e docs: developer pages (integration, mass lumping), DIVERGENCES (C82 rationale, TRIA6-HRZ row), SEGUIMIENTO register, BACKLOG C85 closed.
- [x] T6f C84 `post_element_force` 2D/3D on non-tensor elements and the 2D record (BACKLOG C84). Done `f2ae525` (fix + 14 tests + build_safe), `651de63` (docs). Route: delegated writer (trigger: parent long-session context backstop; calcul_force.cc + tests + build_safe + docs, Pro probes, omen corpus A/B). Risk: item 3 only if an existing input uses the 7-value 2D record (search before the change); the 2D tensor path changes (NaN today), the 3D hex8/hex27 path must stay byte-identical. Subtasks:
  - [x] C84a manual 6.927 record/result layout; corpus search for 7-value 2D records; Pro probes per element type (tria3, tria6, quad9 2D; tet4, tet10, prism6, prism15, hex27 3D).
  - [x] C84b RED tests under test-2026 (tet4 split of `tpef_force10`, 2D quad4 and tria3 strip with 6 values, wrong-count expected error; Pro targets 4e-6) + build_safe registration.
  - [x] C84c implementation: 2D record 6 values (wrong count = input error), 2D stress layout, non-tensor elements through pol() with the real point count; unsupported types refused.
  - [x] C84d verification: clean build, build_safe TN_CHECK_ONLY, audit_db_buffers, hex8 3D byte-identical, omen corpus A/B.
  - [x] C84e docs: user/dev post_element_force (+ README), DIVERGENCES if a row applies, SEGUIMIENTO register, BACKLOG C84 closed.
- [x] T6f LOW queue closure (owner 2026-10-17: fix the cheap items of the low-severity queue, document the rest as class B deferred, then merge and move to F5). Route: delegated writer (trigger: parent long-session context backstop). One work-unit commit per unit with its tests and docs. Subtasks:
  - [x] Q1 PEF 1D result length (Pro `force16`/`post7` record length; calcul_force.cc 1D branches, database.cc). Route: inline in the writer. Risk: item 3 (result record length of an existing output; no reader of index 1-4 in the code).
  - [x] Q2 C60 STOPPED before any code change: not output-only (the GNU orientation rule makes two of four node orders treat a compressed MC/tension interface as open); re-scoped in BACKLOG C60, owner decision.
  - [x] Q3 C62 `parallel_new_dof_before` in the post-solve pass (top.cc), double-application check, dedicated omen A/B of `node_rhside`; re-measure C63. Risk: item 3 (reported `node_rhside` output) + item 6. STOP if any case moves away from Pro.
  - [x] Q4 C54 `db_type()` check of run-time record names read into `idum[1]` (post.cc, print_el.cc, repeat.cc, data.cc, calcul.cc), expected-error test, audit_db_buffers. Risk: item 3 (new input error for wrong-typed names).
  - [x] Q5 C84-REFUSE STOPPED: the Professional INCLUDES a truss with a node on the section (analytically correct total, 100 = 45 + 55), neither ignored nor wrong; new BACKLOG C86, owner decision. No code change.
  - [x] Q6 C58 element history records after `control_restart` (interface, contactspring, spring, tendon, damage) vs Pro, one probe per family; `number_timesteps` -> class B; time box half a day (interface + spring first). Risk: item 1 (restart state).
  - [x] Q7 bookkeeping docs: class B deferrals (C61, C65, C69, C68b, C68c, K14 hypo3 hyhis residual, C22, C12, C25), stale lines (C3/C5 headings, C51 row, DIVERGENCES K21, AUDIT K14/K15/K21, task lines T6/T6f), BACKLOG open recount. C68a / C72-V1 class D and sfnet B14 untouched (owner decisions).

## Later (owner 2026-10-10, after convergence with ALL divergences documented and accepted)
- F6: ILU/AMG preconditioning of the iterative solver (BACKLOG C81, owner 2026-10-15: in large meshes with many steps even small per-solve gains compound).
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
  T6f C74 (2026-10-13), route: delegated writer (trigger: Pro/GNU probes,
  corpus A/B on omen, so.cc + so_bicg.cc + headers, test, 5 docs).
  Investigation: `so_bicg.cc` Bi-CG breakdown test `|r1.r2|`,
  `|d2^T A d1| < 1e-16` is absolute (the CG one relative since
  2026-09-13); near the validation_1 steady state it fires on healthy
  recurrences (normalised r1.r2 1e-3..1e-4, dAd 1e-2..1e-4) in 599 of
  600 solves and `so.cc` re-solves exactly with the band LU (= SuperLU).
  38 of 956 corpus runs printed the error block. Pro 02-08-2026: default
  PARDISO, its Bi-CG prints nothing. A scale-relative variant removes all
  breakdowns but loses accuracy on the absolute 1e-16 floor (examp16
  node_dof vs SuperLU 0.43 % -> 2.45 %, 1e-5 scaled validation_1 8 %; Pro
  Bi-CG the same floor behaviour); a restart variant failed. Owner
  decision 2026-10-13: option A now, B -> C77. Fix `cf07b19`: the block
  moves to `solve_iterative_bicg_breakdown_report()`, printed only when
  the LU retry fails (or on the displacement-primary second solve); a
  rescued breakdown prints "Warning: Bi-CG breakdown, system solved with
  the direct band LU". Test test-2026 bicg_breakdown_lu_retry1 + build_safe
  check: RED on 46104a3 (rc 0, 9 error blocks, check FALLO), GREEN (9
  warnings, no block). Clean build; build_safe TN_CHECK_ONLY=1 263 OK /
  0 FALLO, runs 226/226; audit_db_buffers exit 0. Corpus A/B omen
  c74b_runs (957 runs): no classification change, 881/882 non-timeout
  .dbs identical masking time_calculation/time_at_start (the other,
  bicg_csr_threads1, 2-thread round-off that also differs between HEAD
  runs); error-block runs 38 -> 0. audit_all_variables.py FATAL regex
  unchanged (message text unchanged). Docs `e5f50fc`: BACKLOG C74
  resolved + C77; DIVERGENCES C77 class B (B 45, total 89); SEGUIMIENTO
  register; dev solve_iterative_bicg (stale 1e-10/1e-12 defaults fixed to
  1e-14/1e-16, breakdown + LU retry), user solver. Next: C57 post_point,
  C71, K21 validation_1, C76, C77.
  T6f C57 (2026-10-13), route: delegated writer (trigger: Pro/GNU
  probes, corpus A/B on omen, mesh.cc + test + 9 docs). Cause:
  `mesh.cc` `mesh_extrude()` quad9 -> hex27 branch (C47 shared planes)
  copied `node_start_refined` from the z = 0 source node, so the inner
  plane nodes kept z = 0 there and `point_el()` (which reads
  `NODE_START_REFINED`) located post points in collapsed hex27. Pro
  02-08-2026 post_point = exact hex27 interpolation of the nodal dofs
  (1e-14). Other mesh branches (mirror, copy, rotate 2d/3d, move,
  switch, convert_quad8/hex20, linear extrude) checked: already correct.
  Fix `cf37de3`: rewrite `NODE_START_REFINED` with the plane z. Test
  test-2026 post_point_extrude_hex27_1 (4 post points, Pro targets, abs
  1e-11): RED on 8b36a34 rc 1 (disy -3.51050e-5 / -1.13556e-5 /
  -3.25367e-6 / -1.78447e-5), GREEN rc 0 (point 1 9.6e-12 off: the
  point_el descent stops ~1e-8 iso short). Clean build; build_safe
  TN_CHECK_ONLY=1 264 OK / 0 FALLO, runs 226/226; audit_db_buffers exit
  0. Corpus A/B omen c57_runs (958 runs): no classification change,
  878/882 .dbs identical; force11 / extrude_quad9_shared_nodes1 differ
  only in node_start_refined (7200 / 45 nodes, now = node), plus the new
  test and bicg_csr_threads1 round-off. Curved quad9 / hex27 probes: Pro
  locates with the corner nodes only (bilinear / trilinear), GNU exact
  isoparametric; owner 2026-10-13: deliberate divergence C78 class A.
  Docs `30f2621` (C57 cause, C77 owner notes), `1f5e82c`: BACKLOG C57
  resolved, C78 decided (49 open); DIVERGENCES C78 class A (A 27, B 44,
  total 89), C57 removed; SEGUIMIENTO register + Descarte rows;
  user/developer post_point (new) and control_mesh_extrude. Next: C71,
  K21 validation_1, C76, C77.
  T6f C71 (2026-10-13), route: delegated writer (trigger: Pro/GNU
  probes, corpus A/B on omen, polynom.cc + test + 6 docs). First NaN:
  `polynom.cc` `pol()` axisymmetric hoop entry `new_b += h/radius` at the
  Lobatto points on r=0, times the zero point volume 2*PI*r*w*detj ->
  NaN velx node_rhside -> NaN velocities. Pro 02-08-2026 probes
  (element_intpnt_dof): an element touching the axis in a -lobatto
  axisymmetric group = its -gauss run (quad4 and quad9), off-axis
  elements keep Lobatto; switch absolute on the radius between 9e-10 and
  1e-9 (EPS_COORD). Manual 6.621/6.623 silent. Fix `bea7599`: pol()
  rebuilds the rule with Gauss (and stress_recovery_gauss=1) when an
  axisymmetric Lobatto point has |r| < EPS_COORD. Test test-2026
  axisym_lobatto_axis1 (Pro targets 1e-9): RED rc 1 on 58af21e (node 3
  vely 1.1 vs -0.96, node 1 stresses 0), GREEN rc 0. sfnet tube1
  (modern keywords): velocities = Pro within 7e-9; nodal stress at the
  r=1.25 column off (-494.42 vs -407.64): Lobatto quad4 recovery is
  volume weighted in GNU, arithmetic mean in Pro (two plane quad4 of
  widths 1/2: 6.525 vs 6.15) -> new BACKLOG C79 (not fixed: other
  Lobatto rules unmeasured). Clean build; build_safe TN_CHECK_ONLY=1 265
  OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Corpus A/B omen
  c71_runs (959 runs): only tube1 (closer to Pro 0.59 -> 0.18) and the
  new test changed; ground11_nonsaturated 118.7 s -> timeout edge (Pro
  also times out). Pending: calcul_force.cc msf_element_rule() replica
  (outside the surface). Docs `db627e8`: BACKLOG C71 resolved + C79 (49
  open); DIVERGENCES C71 row -> C79 (B 44, total 89); SEGUIMIENTO
  register; user/dev integration_quadrature_and_stress_recovery + README
  indexes. Next: K21 validation_1, C76, C77, C79.
  T6f C79 (2026-10-13), route: delegated writer (trigger: Pro/GNU
  probes, corpus A/B on omen, miscel.cc + polynom.cc + calcul_force.cc +
  headers + 6 tests + 9 docs). Pro 02-08-2026 probes (12 probes,
  prescribed nonlinear velocity at every node, distorted meshes of
  unequal element size): Lobatto with the points at the nodes (quad4 2x2,
  quad9 3x3, hex8 2x2x2, hex27 3x3x3; plane + axisymmetric off-axis) =
  point value at the node + arithmetic mean (5e-12; physical fit
  0.39-0.46 off); quad9 2x2 / hex27 2x2x2 Lobatto -minimal = physical
  linear fit of the corner points + mean (5e-12; point value 0.18-0.19
  off, also on a rectangle); quad4/hex8 1-point already = Pro. Manual
  6.623 silent -> Pro oracle, no product decision. GNU cause:
  sri_stress_recovery_weight() returned h for every Lobatto rule (x point
  volume = volume weighting; never measured, 06e175d kept h); quad9/hex27
  -minimal non-corner nodes got nothing (materi.cc h guard; sigxx 0 vs
  up to 16.5). No earlier test/corpus target used -lobatto. Fix
  `ab4ef78`: recovery_fit_kind() kind 3 (point value) for npoint==nnol
  Lobatto quad4/quad9/hex8/hex27; kind 1 fit with Lobatto point
  positions for the -minimal ones; pol() sets stress_recovery_gauss=1 for
  Lobatto with fewer points than nodes and the new
  stress_recovery_lobatto (both thread_local). C71 leftover:
  msf_element_rule() rebuilds the rule with Gauss when an axisymmetric
  Lobatto point has |r| < EPS_COORD. Tests test-2026
  recovery_lobatto_{quad4_mean1, quad4_axisym_mean1, hex8_mean1,
  quad9_minimal1, hex27_minimal1} (Pro targets 4e-6) and
  msf_axisym_lobatto_axis1 (target = -gauss run momx_sig 1/6; Pro prints
  zero section forces for this axisymmetric input): RED rc 1 on d2a365f,
  GREEN; all probes every node = Pro within 5e-11; tube1 every node_dof =
  Pro within 6e-8. Clean build; build_safe TN_CHECK_ONLY=1 271 OK / 0
  FALLO, runs 226/226; audit_db_buffers exit 0. Corpus A/B omen c79_runs
  vs c71_runs (959 runs): tube1 node_dof -> Pro; slide2, slide3,
  axisym_lobatto_axis1 only stress-dof node_rhside (now = Pro, node_dof
  unchanged); bicg_csr_threads1 (not an identity witness). Interface
  elements not involved (interface.cc assembles forces only). Found C80:
  group_integration_method -reduced refused by the GNU (Pro accepts it).
  Docs `6094739`: BACKLOG C79 resolved + C80 (49 open); DIVERGENCES C79
  row -> C80 (B 44, total 89); SEGUIMIENTO register; user/dev
  integration_quadrature_and_stress_recovery, nodal_stress_recovery,
  post_calcul_materi_stress_force. Next: K21 validation_1, C76, C77, C80.
- 2026-10-15 T6f C77 (route: delegated writer; trigger: so_bicg.cc +
  so.cc + test + build_safe + 5 docs, corpus A/B on omen, local timing).
  Fix `d567672`: relative stop |r|^2 < bicg_error |r0|^2 (CG and Bi-CG),
  bicg_error_minimum only from the user record (GNU record; not in the
  Pro or 2014 manuals), round-off floor 100 eps |A||x|, relative Bi-CG
  breakdown 1e-14, floor/breakdown iterate kept at |b-Ax|/|b| <= 1e-6,
  else SuperLU retry then band LU within 0.7 x free memory (clear error,
  no bad_alloc); no stagnation window / best iterate. Test
  bicg_relative_stop1 (validation_1 inflow x1e-5): RED HEAD rc 1 (4.0e-5
  off SuperLU), GREEN 2e-10. bicg_breakdown_lu_retry1 check: no false
  breakdown (target unchanged). Clean build; build_safe 272 OK / 0 FALLO,
  runs 226/226; audit exit 0. 86 Bi-CG inputs: retries 7378 band -> 13
  SuperLU, 19 closer / 0 farther, wall +3.5 %. Large: C23 91.9 -> 18.6 s,
  C36 bad_alloc -> 83 s, CG E25 +7 % wall, 9.3e-7 -> 1.6e-11. Omen A/B
  c77_runs/rs vs c79_runs (965): ground11_nonsaturated timeout -> rc 1
  (its rc), truss13 rc 1 both (singular, now a reported breakdown); 84
  moved runs: 63 toward SuperLU (c77_runs/*/slu), 10 same, 1 away
  (dynamic6 diverging), 10 other refined mesh. BiCGStab: large -9..13 %
  solver time, corpus 71 retries / 29 max_iter / 2 farther; GMRES(30)
  clearly worse -> question to the owner (Bi-CG kept). Docs `5091fe7`:
  DIVERGENCES C77 class A (A 28, B 43), SEGUIMIENTO register + Descarte
  MKL PARDISO, BACKLOG C77 resolved + C81 (F6). Next: owner BiCGStab
  choice; K21 validation_1, C76, C80.
- 2026-10-16 T6f C80 (route: delegated writer; trigger: headers +
  database.cc + polynom.cc + miscel.cc + calcul_force.cc + delete.cc + 7
  tests + build_safe + 9 docs, corpus A/B on omen). Owner decisions
  2026-10-16: default factor 10 (manual 6.621; Pro applies 1) class A;
  hourglass warning in the user manual; selective/B-bar -> C80.1 (OPEN,
  next, not implemented). Fix `995de49`: enums REDUCED +
  GROUP_INTEGRATION_METHOD_REDUCED_FACTOR (both headers), database.cc
  name + record (DOUBLE, length 1, required GROUP_TYPE), polynom.cc
  integration_reduced_factor() (default 10, f <= 0 db_error) and pol()
  quad4/hex8 points +-1/f weight 1/2, other cases Gauss/Lobatto fallback,
  quad9/hex27 Lobatto -> kind 1 fit, stress_recovery_reduced_iso for the
  miscel.cc fit points; axisymmetric node on the axis -> Gauss (new Pro
  probe: f=1 AND f=2 switch on the NODE radius); msf_element_rule /
  mc_element_rule replicas. Tests test-2026 reduced_quad4_f2_1,
  reduced_quad4_f10_1, reduced_hex8_f2_1, reduced_quad9_minimal1,
  reduced_quad9_maximal1, reduced_quad4_axisym1 (Pro 4e-6) and
  reduced_default_factor1 (GNU default = Pro f=10 run): RED rc 1 on
  0ac600d ("I do not know : -reduced"), GREEN; every node_dof = Pro
  within 5e-12; earlier 135-probe set = Pro except bar3 recovery
  (1.3-2.3 %, C67) and prism6 -minimal (C76). Clean build; build_safe
  TN_CHECK_ONLY 279 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0.
  Omen A/B c80_runs/fix vs c77_runs/rs (966 runs): 0 changed (masked
  incl. tochnog_version), only-new the 7 tests rc 0. Docs `6c8269f`:
  user/dev group_integration_method_reduced (new), integration pages,
  README indexes, DIVERGENCES C80 class A (A 29, B 42), SEGUIMIENTO
  register + Descarte, BACKLOG C80 resolved + C80.1 + C67 note. Not
  investigated: Pro element-average fallback at f=1e6. Next: C80.1
  (keyword decision by the owner), K21 validation_1, C76, C77c.
- 2026-10-16 T6f C80.1 (route: delegated writer; trigger: miscel.cc +
  materi.cc + database.cc + calcul_force.cc + headers + 4 tests +
  build_safe + 6 docs, omen corpus A/B). Owner decision 2026-10-16:
  extend `group_element_selective_reduced_integration` with -volumetric
  / -both (GNU-only, opt-in). Formulation: D = D_dev + kappa m m^T,
  kappa = m^T D m / 9, volumetric part at the centroid point, D-c
  feedback + B-bar stored stress; sri_active() mode mask (volumetric
  dropped on membrane groups, invalid value db_error); calcul_force
  shear replacement only for the shear bit; check.cc not touched (no
  value checks there). Fix `b47c39c`. RED rc 1 on 66df7e1 ("I do not
  know : -volumetric"/"-both"), GREEN. Lame 8x16 nu 0.4999 -0.19 %
  (Gauss -80 %, Pro f=10 -11 %), checkerboard 4 % vs 12531 %; hex8 =
  quad4; no hourglass; quad4 3 / hex8 6 zero eigenvalues (-both hex8 9).
  Coarse bending soft (cantilever 10x1 nu 0.3: -volumetric +22 %, -both
  +113 %; exact element behaviour, an independent displacement solve
  gives the same digits); lambda/mu split variant measured better in
  coarse bending (-4.9 %) but not chosen (owner question). Gauss/-yes/
  default study runs byte-identical. Clean build; build_safe
  TN_CHECK_ONLY 283 OK / 0 FALLO, runs 226/226; audit_db_buffers exit
  0. Omen A/B c801_runs vs c80_runs (973 runs): 0 node_dof changes;
  bicg_csr_threads1 post_point 1e-15 text change = threaded round-off
  (baseline varies run to run); only-new the 4 tests rc 0. Docs
  `7504e47`: user/dev SRI pages, READMEs, SEGUIMIENTO, BACKLOG C80.1
  resolved (48 open); DIVERGENCES unchanged (shear SRI not listed).
  Open owner questions: lambda/mu split, B-bar for plasticity,
  axisymmetry, adaptive reduced factor.
- 2026-10-16 T6f C80.2 (route: delegated writer; trigger: miscel.cc +
  materi.cc + database.cc + input.cc + headers + 4 tests + build_safe +
  docs, omen corpus A/B). Owner decision 2026-10-16: `-lambda` (Hughes
  1977 lambda/mu split, isotropic linear elastic), `-both` removed, gate
  before committing. Implementation: C80.1 machinery with kappa ->
  lambda (SRI_LAMBDA 4), lambda = mean off-diagonal of the normal block,
  isotropy check (1e-6) else error naming -volumetric; plane stress
  dropped; enum BOTH -> LAMBDA in place; input.cc dedicated -both
  error. Gate PASSED: Lame 8x16 nu 0.4999 -0.24 % (vol -0.19 %),
  sigma_theta(a) -0.14 % (vol -3.1 %); element pressure oscillation
  0.1 % (= vol; IP intra-element 13 %, nu-independent, O(h)); TG
  cantilever / pure bending 10x1 nu 0.3 -4.9 % (vol +22 %), nu 0.49
  +30 % (vol +67 %), 40x4 nu 0.49 +1.4 % (vol +2.5 %); hex8 = quad4;
  distorted 20x2 lambda -12/-3.8/-3.3 % vs vol -8.4/-0.7/-0.3 % (no
  locking trend, documented); hourglass 1.5x Gauss, zig-zag 0.083;
  zero eigenvalues quad4 3 / hex8 6 also distorted. Independent solves
  match GNU (2D distorted 1e-10, 3D distorted 1e-12). RED on 051c706
  (rc 1 "I do not know : -lambda"; -both rc 0), GREEN. Clean build;
  build_safe TN_CHECK_ONLY 286 OK / 0 FALLO, runs 226/226;
  audit_db_buffers exit 0. Omen A/B c802_runs vs c801_runs (977 runs):
  0 node_dof changes; text only bicg_csr_threads1, alternate1 (timeout
  border), cylinder_buckle (timeout both); sri_both_cantilever1 now the
  -both error. Docs `7866a53`: user/dev SRI pages, READMEs, SEGUIMIENTO,
  BACKLOG (C80.2 note on C80.1). Open: build_safe expected_error_ok
  greps the marker without `--`.
- 2026-10-16 T6f C67 (route: delegated writer; trigger: miscel.cc +
  polynom.cc + 10 tests + build_safe + 9 docs, omen corpus A/B).
  Professional measured on 65 probes (`scratchpad/c67/gen67.py`,
  `an67.py`; bar2 x4 / bar3 x3 / bar4 x3 elements in 1D, quad16 2x2;
  default/-gauss/-lobatto/-reduced f=2,4 x default/-minimal/-maximal),
  every node: Gauss, -reduced and Lobatto -minimal = physical linear fit
  1,x(,y) (5e-13; natural fit 1e-2, quadratic 5e-2); Lobatto -maximal =
  point value (bar4/quad16 interior nodes at the Lobatto positions in
  both codes, Pro element_intpnt_h = identity); arithmetic mean
  (length-weighted 1e-2 off); quad16 axisymmetric on the axis = C71 Gauss
  switch + sum h sig / sum h. hex64 rejected by both codes (Pro "I do
  not know : -hex64", GNU element length 28 < 65); truss/beam have no
  recovered nodal stress. GNU cause: h weights (unmeasured fallback, not
  intended; bar3 Lobatto -minimal mid nodes 0); no earlier test covered
  these rules. Fix `d875874`: recovery_tensor_npol(), dimension-generic
  recovery_fit_kind()/recovery_fit_shape() (nodes at Lobatto positions),
  pol() fits every -reduced Lobatto rule; point rules unchanged, so no
  msf_element_rule/mc_element_rule change. Tests test-2026
  recovery_{bar2_gauss2_mean1, bar3_gauss_fit1, bar3_lobatto_minimal1,
  bar3_reduced1, bar4_gauss_fit1, bar4_lobatto1, quad16_gauss_fit1,
  quad16_lobatto1, quad16_lobatto_minimal1, quad16_axisym_lobatto1}: RED
  rc 1 on 75bc0b5, GREEN; 56 supported probes every node = Pro within
  4e-12. Clean build; build_safe TN_CHECK_ONLY 296 OK / 0 FALLO, runs
  226/226; audit_db_buffers exit 0. Omen A/B c67_runs vs c802_runs (981
  runs): 0 node_dof changes; text only alternate1 (timeout border),
  cylinder_buckle (timeout both), bicg_csr_threads1; only-new the 10
  tests rc 0. Docs `aee381e`: user/dev nodal_stress_recovery,
  integration_quadrature_and_stress_recovery,
  group_integration_method_reduced; DIVERGENCES (closed-by-following
  list, C80 row), SEGUIMIENTO register, BACKLOG C67 narrowed to tet4
  -maximal, prism15, simplex axisymmetric (48 open).
- 2026-10-16 T6f C67 remainder + C76 (route: delegated writer; trigger:
  Pro/GNU probes + miscel.cc + polynom.cc + 5 tests + build_safe + 6
  docs, omen corpus A/B). Enumeration of every `pol()` rule against the
  recovery dispatch: on the `h` fallback were tet4 `-maximal`, prism15
  (all values) and prism6 `-minimal` (record ignored). 81 Professional
  probes (`scratchpad/c76/gen76.py`, `an76.py`; tria3, tria6, tria3
  axisymmetric off/on axis, tria6 axisymmetric, tet4, tet10, prism6,
  prism15 x 3 methods x 3 point values). Pro: tet4 `-maximal` = Keast 5
  (manual 6.623) + linear fit; prism15 = linear fit (5e-13); prism6
  `-minimal` = 2 centroid points, `sum h sig / sum h` (5e-13); all
  arithmetic mean. Found a GNU prism15 rule bug: triangle points
  (b,b)/(d,d) instead of (a,a)/(c,c) (one outside the triangle,
  monomial error 0.025). C76 analytic: 2-point rule 1 spurious mode per
  element, restrained by every mesh tried, patch test exact, bending
  -7.1 % vs -13.1 % (20x2x2) -> sound, implemented; weights confirmed
  (cantilever = 2-point FE model to 1.5e-13). Pro prism15 `-minimal` (8
  points, centroid -27/48) has 3 negative stiffness eigenvalues -> not
  implemented, C82 (owner). tria6: Pro = manual 7/4 points + fit, GNU 6
  nodal (C66 re-scoped: lumping decision in general.cc needed). tria3
  `-maximal`: manual 4, both 3; tet10 `-maximal`: manual 10, both 5;
  axisymmetric tria3 `-maximal` GNU 1 point vs Pro 3 nodal / on-axis 1
  -> C83. Fix `dc4a59a`: pol() tet4 Keast, prism6 2-point (default set
  to -maximal first: a regression caught on the probes), prism15 table;
  recovery_fit_kind/shape/build (kind 4); replicas unchanged (tensor
  only). RED rc 1 at 9a0ca3d (5 tests), GREEN; every probe node of a
  changed rule = Pro within 6e-12. Clean build; build_safe TN_CHECK_ONLY
  301 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Omen A/B
  c76_runs vs c67_runs (991 runs): only `prism15` changed (node_dof vs
  Pro 2.1e-14 -> 3.3e-15, toward Pro), `bicg_csr_threads1` text only
  (nondeterministic, node_dof 0); only-new the 5 tests. Docs `5bd6f90`:
  user/dev nodal_stress_recovery and integration pages, prism15 page;
  DIVERGENCES C76 closed, C66/C82/C83 class B (B 44, total 91);
  SEGUIMIENTO register + group_integration_points line; BACKLOG C67/C76
  resolved, C66 re-scoped, C82/C83 new (48 open).
- 2026-10-16 T6f C82 + C83 (route: delegated writer; trigger: probes +
  polynom.cc + test + build_safe + 9 docs, omen corpus A/B). Owner
  decisions: C82 and C83 part 1 class A (docs only), C83 part 2 follow
  the Professional. Why the GNU forced 1 point: `pol()` line from the
  2016 initial import, no comment; the only reason found is the r=0 hoop
  singularity (C71), avoided by the per-element on-axis switch. Pro
  (C76 probes `scratchpad/c76/`): axisymmetric tria3 `-maximal` = 3
  nodal points, 1 centroid point in the elements with a node on the axis
  (elements 1-4 of the ax0 probe 1 point, 5-12 3 points); default and
  `-minimal` 1 point (= GNU 3e-12). Fix `72f78bd`: pol() tria3 on-axis
  switch (EPS_COORD on the node radius); recovery unchanged; replicas
  tensor-only (not touched). RED rc 1 on 55e6a10, GREEN; 27 tria3 probes
  = Pro within 3.1e-12 (before 3.2e-2). Clean build; build_safe
  TN_CHECK_ONLY 302 OK / 0 FALLO, runs 226/226; audit_db_buffers exit
  0. Omen A/B c83_runs vs c76_runs (996 runs): 0 changed (no corpus case
  has an axisymmetric tria3 `-maximal`); only-new the test rc 0. tet10:
  a positive 10-point degree-3 rule exists (Keast: 4 x (a,b,b,b)
  a=0.5684305842 w=0.2177650699 + 6 edge midpoints w=0.0214899534),
  not implemented. Docs `8e9606d`: DIVERGENCES C82/C83 class A (A 31,
  B 42, total 91), SEGUIMIENTO Descarte x2 + register, user/dev
  nodal_stress_recovery, integration pages, prism15; BACKLOG C82/C83
  closed (46 open). Open: `post_element_force` (`calcul_force.cc`) feeds
  a tria3 to the quad4 msf helpers (pre-existing, unrelated).
- 2026-10-16 T6f C66 phase 1 (route: delegated writer; trigger: Pro/GNU
  mass probes + lumping-site inventory). Owner: measure the Professional
  lumping first. Method: one element, density 1, E 1e-9, unit point
  force at one node, every other node fixed, one step dt 0.01, m = F dt
  / v (`scratchpad/c66/genm.py`, `anm.py`; 174 Pro + 174 GNU runs);
  condif capacity cross-check (`scratchpad/c66/cd/`). Pro = row sum of
  the consistent mass for every element and rule value: tria3 A/3, quad4
  A/4, quad9 A/36 / A/9 / 4A/9, tet10 -V/20 / V/5 (GNU = Pro except
  tet10 HRZ); tria6 plane corner 0 (free one-element step: Pardiso zero
  pivot), mid-side A/3; tria6 axisymmetric (r 1..2) corners -pi/60,
  -pi/60, +pi/30, mid-side 1.466/1.257/1.466 (total 2 pi A r_c); capacity
  identical. GNU HEAD tria6 corner 5e-7. Gate: zero/negative corner
  masses -> stopped before any source edit (C66d not started).
  Inventory: general.cc inertia block h-weight (principal dofs + the
  h-weighted state dofs) affected; recovered dofs (F/volume fit weights)
  not affected. Docs `9b51253`: BACKLOG C66 measured lumping + owner
  question (HRZ for tria6 vs keep the GNU rule), new C84
  (`post_element_force` 2D/3D: tet4 cantilever normal 3143 vs Pro
  -123.4; 2D record 7 vs 6 values and NaN), 47 open; DIVERGENCES C66
  row. Next: owner decision on C66.
- 2026-10-16 T6f C66 phase 2 (route: delegated writer; trigger:
  polynom.cc + miscel.cc + general.cc + materi.cc + headers + 5 tests +
  build_safe + docs, Pro probes, omen A/B). Pro: tria6 iso = L1,L2,
  7-point Hammer (10-digit literals) / 4-point Strang-Fix; on-axis
  axisymmetric tria6 (new probes `scratchpad/c66/p2/probes/tria6ax0_*`):
  no rule switch, linear fit (no h-average), no NaN. Fix `175faf1`:
  pol() 7/4-point rules + `tria6_shape()`, axisymmetric refusal removed,
  thread-local `lumped_hrz_weight` (7-point consistent diagonal scaled
  to the total, constant per point, rule-independent masses); general()
  HRZ for every h-weighted tria6 dof except groundflow_velocity
  (groundfl.cc numerator out of surface -> BACKLOG C85, with condif
  absorption); materi() `h_lump` state numerators; recovery_fit_kind/
  shape kind 1, no on-axis h-average for tria6; replicas tensor-only.
  RED rc 1 on c730f39 (5 tests), GREEN; probes = Pro within 3.1e-10;
  tria6_nodal_recovery_mean1 tightened (now 5e-10). Analytic (bar modes
  1-5 vs (2n-1)pi c/2L): consistent +0..+0.19 %, HRZ -0.01..-0.67 %,
  row sum condensed -0..-0.18 %; time-domain tip peaks identical. 4-point
  stiffness on curved elements: no spurious modes. Clean build;
  build_safe TN_CHECK_ONLY 307 OK / 0 FALLO, runs 226/226;
  audit_db_buffers 0. Omen A/B c66_runs vs c83_runs (997): changed
  `tria6` x2 (3.3e-10, = Pro as before), `tria6_nodal_recovery_mean1`
  (toward Pro); slope_classical_numerical (timeout both, rerun 1800 s):
  collapse 1.26297 -> 1.26265 vs Pro 1.27045, final stresses vs Pro 20.8
  -> 0.25. Docs `4b65b07`: DIVERGENCES TRIA6-HRZ (A 32, B 41, total 91),
  SEGUIMIENTO register + Descarte, BACKLOG C66 closed + C85 (47 open).
  Next: C85 (needs groundfl.cc / condif.cc in a surface).
- 2026-10-17 T6f C85 + C82 consistency (route: delegated writer;
  trigger: groundfl.cc + condif.cc + general.cc + 3 tests + build_safe +
  docs, eigenvalue study, omen A/B). Fix `50fd4c7`: `h_lump` in the
  groundfl.cc velocity numerator and the condif.cc absorption (rhs,
  lhside, matrix, residue); general.cc groundflow_velocity exception
  dropped. Tests (hand targets) RED rc 1 on c6fa9aa, GREEN:
  `tria6_groundflow_velocity_hrz1` (gvelx 10 everywhere; before corners
  35.2/-40/-144/-368), `tria6_condif_absorption_hrz1` / `_axisym1` (1/1.1
  everywhere; before corners 0.9955/1.0027, axisym up to 1.137). Clean
  build; build_safe 310 OK / 0 FALLO, runs 226/226; audit_db_buffers 0;
  omen A/B c85_runs vs c66_runs (999 common): changed only
  bicg_csr_threads1 (nondeterministic), only-new the 3 tests. Docs
  `eb00f27`: dev integration + inertia_apply, user integration,
  DIVERGENCES TRIA6-HRZ row, SEGUIMIENTO register, BACKLOG C85 closed
  (46 open). C82 check (`scratchpad/c85/eig_tria6.py`, `eig_extra.py`):
  tria6 4-point (-27/48) stiffness exact on straight plane elements
  (|K4-K7| 7e-16; 3 zero, 0 negative at nu 0.3/0.49); ONE negative
  eigenvalue on curved plane elements at nu 0.49 (min/max -4.4e-3 /
  -2.9e-2; onset nu 0.487 / 0.434) and axisymmetric at nu 0.49 (straight
  on-axis -1.3e-2, onset 0.386; straight r 1..2 onset 0.496; curved r+1
  -4.2e-3; r+10 none); 7 and 3-point rules PSD everywhere. prism15
  8-point: 1 negative + 2 spurious zeros already on a regular prism at
  nu 0.3 (C82 recorded 3 negatives on its distorted probe). Stopped,
  owner question: keep the tria6 4-point rule (= manual / Pro) or align
  with C82. Next: owner decision on tria6 `-minimal`.
- 2026-10-17 owner decision (tria6 `-minimal`): keep the manual/Professional
  4-point rule; its indefinite stiffness on curved/axisymmetric elements with
  nu > ~0.4 is documented as a user-manual limitation (integration page);
  criterion vs C82: refuse a rule that fails in normal use (prism15, regular
  prism at nu 0.3), keep with a warning one that fails only in a corner case.
  C77c ticked (owner 2026-10-15: keep Bi-CG). SEGUIMIENTO global metric
  refreshed (1154 / ~833 accepted / ~258 pending). Next: C84
  (`post_element_force` element types and 2D record), then the LOW-queue
  triage, then merge + F5.
- 2026-10-17 T6f C84 (route: delegated writer; trigger: parent
  long-session context backstop). Manual 6.927: 2D record 6 values
  (n_x n_y s0_x s0_y mid_x mid_y), 2D result normal / shear / moment;
  no corpus, test-2014, test-2026, sfnet or Pro `other` input uses a
  7-value 2D record (Pro force14/force15 use 6), so the record change
  is not a contract break. Pro 02-08-2026 probes (`scratchpad/c84/`):
  tpef_force10 split into tet4/tet10/prism6/prism15 -> -123.4 / 0 /
  1000 / 0 / -100000 (GNU HEAD: normal 3143.7 / -1651.5 / 4918.6 /
  1871.5); 2D strip quad4/tria3/quad9/tria6 -> 3 values 100 / 10 / -40
  (GNU: 6 values skipped, 7 values NaN: the 2D stress copy wrote sxx
  sxy sxz into a 6-stride array read with stride 3 as sxx syy sxy);
  axisymmetric ring strip -> normal 100 (ring forces, hoop term); 7 or
  5 values -> "Number of data values expected : 6". Fix `f2ae525`:
  input.cc count check, pol() path `pef_element_internal_forces_pol`
  for the non-tensor elements and axisymmetric 2D, 2D layout fix,
  3-value 2D result (fixed_length 0), refusal of non-continuum /
  interface elements; found in the omen A/B: interior sections summed
  both sides (cancel; force14 shear 33.3 vs Pro -500) -> mean of the
  positive side and the negated negative side (force14/force15 =
  Pro exactly, reversed normal = Pro). Tests 14 (RED on 54b8170, GREEN);
  hex8/hex27 boundary sections byte-identical. Clean build; build_safe
  324 OK / 0 FALLO, runs 226/226; audit_db_buffers 0; omen A/B
  c84b_runs vs c85_runs (1005): force14/force15 rc 1 -> 0 (= Pro),
  bicg_csr_threads1 (nondeterministic), only-new the 14 tests. Docs
  `651de63`: user page (new) + dev page, DIVERGENCES class B
  C84-REFUSE (B 42, total 92), SEGUIMIENTO register, BACKLOG C84
  closed (45 open). Open: Pro handling of trusses/interfaces in the
  selection (C84-REFUSE), 1D result length. Next: LOW-queue triage.
- 2026-10-17 T6f LOW queue closure (route: delegated writer; trigger:
  parent long-session context backstop; owner plan 2026-10-17: fix the
  cheap low-severity items, defer the rest as class B, then merge + F5).
  Q1 `feafe50`: 1D `post_element_force_result` = one value (Pro
  force16 5.034285714286; target on value 1 -> target_item error); tests
  `post_element_force_1d_length1` / `_length_error1` (RED on b938b0f).
  Q2 C60 STOPPED: not output-only. Probes `scratchpad/lowq/c60/`: blocks
  + quad4 interface, four node orders: Pro -9.98 always; GNU +/+/-/- and,
  with MC + tension cut-off, orders `4 2 10 9` / `10 9 4 2` treated as
  OPEN (normal 0). Free interface: Pro numbering rule, opened elastic
  interface carries 0. Owner decision (BACKLOG C60 re-scoped).
  Q3 `a2e68b6`: the pass also lost the `bounda_force` loads; top.cc
  restores the iteration's post-bounda() NODE_RHSIDE snapshot and runs
  `parallel_new_dof_before`. Tests node_rhside_post_pass_{node_force,
  node_mass,bounda_force}1 (Pro targets, RED on feafe50). Omen A/B
  lq3_runs vs c84b_runs (1019): node_dof identical everywhere, rc
  unchanged, node_rhside changed in 133 cases; of 64 with a Pro run 62
  only closer to Pro, trubea2/3 (rc 1 both, C37 planar beam, disz 2e6 vs
  Pro 0.033) now show the applied 0.01 at the loaded node (Pro 0): judged
  not a regression (the old 0 hid the load). C63 re-measured: unchanged
  -> class B + DIVERGENCES row.
  Q4 `8edaf3b`: db_type() check in post.cc (post_node, post_integrate),
  print_el.cc, repeat.cc (save, until); data.cc/calcul.cc already safe.
  4 expected-error tests (RED rc 0 / endless repeat). Pro: zeros for
  post_node/post_integrate -> DIVERGENCES C54-TYPE (class C).
  Q5 C84-REFUSE STOPPED: Pro includes a truss with a node on the section
  (100 = 45 + 55, moment -40), correct total -> new BACKLOG C86, owner.
  Q6 `ef9af63`: Pro keeps the interface and spring history over
  control_restart, GNU = Pro; 2 pinning tests; contact-spring, tendon,
  damage, number_timesteps deferred (class B).
  C25 item 1 `efff819` (pef_vel_prev zeroed after reallocation; no
  deterministic RED possible). Q7 `0f707f3`: deferrals (C22, C25 rest,
  C58 rest, C61, C63, C65 -> PPD, C68b, C68c, C69, K14 hypo3 hyhis), C12
  closed (Pro 02-08-2026 does not stop on the singular probe either; GNU
  C77 SuperLU fallback), stale lines (C3/C5, C51, DIVERGENCES K21, AUDIT
  K14/K15/K21, T6/T6f lines), BACKLOG recount: 41 open (2 MEDIUM, 39
  LOW; 32 active, 9 deferred class B; 12 class A decided, 45 closed).
  DIVERGENCES A 32, B 42, C 18, D 1, total 93. Verification: clean
  build; build_safe TN_CHECK_ONLY 335 OK / 0 FALLO, runs 226/226;
  audit_db_buffers 0; omen final A/B lqf_runs vs lq3_runs (1024): only
  cylinder_buckle (timeout both) and the new tests. Open for the owner:
  C60 orientation rule, C86 truss inclusion. Next: merge + F5.
