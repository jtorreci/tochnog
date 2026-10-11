# Backlog — open issues after the convergence cycle

Date: 2026-09-16
Branch: `dev/masin-kernel-surface` @ `e9c40c9`

Every remaining non-PASS corpus test and every residual found during the
GNU→Professional convergence cycle, with its symptom, root cause, evidence,
severity, proposed fix and estimated effort. This is the actionable companion
to `CONVERGENCE-CYCLE-REPORT.md`.

## Severity scale

| level | meaning |
|---|---|
| **CRITICAL** | silently corrupts results on common models, crashes on valid input, or blocks a whole feature family with no workaround |
| **HIGH** | a documented feature diverges from the Professional on a reproducible model a plausible user would hit |
| **MEDIUM** | a narrow or legacy divergence, a coverage gap, or a performance cost that affects a common configuration |
| **LOW** | a single test artifact, a stale/malformed input, or a cosmetic/documentation issue |

## Effort scale

| level | meaning |
|---|---|
| **cheap** | under a day, isolated change |
| **medium** | days, one subsystem, needs an A/B corpus run |
| **large** | weeks, cross-cutting refactor or new infrastructure |

---

## Resolved since the cycle snapshot

These items were open when the cycle snapshot (`a1362ef`) was taken and have
since been closed. They are recorded here so the history is complete; the full
acceptance evidence lives in `CORPUS-STALE-TARGETS.md` and
`SEGUIMIENTO-CONVERGENCIA.md`.

### R1. `iface_tangref` — target recalibrated to the honest per-IP value

- **Corpus**: validation-suite 233.
- **Was**: target `element_interface_force_tang2 2 = 10.4025 ±0.2` (pre-`bf74f28`
  scalar "side average, all-against-all" semantics).
- **Resolution**: since `bf74f28` the record is a length-4 per-integration-point
  array; the honest converged per-IP value is `4.9998` (the interface carries the
  applied load 5; all four IPs agree to 4e-9). The target was recalibrated to
  `4.9998 ±0.2` (tolerance kept), a corpus decision like the 14 earlier
  recalibrations.
- **Commit**: `fbcf144`.
- **Evidence**: `CORPUS-STALE-TARGETS.md` §1 (validation-suite 234: 229/5 →
  231/3; the 233 real tests 229/4 → 231/2; sfnet byte-identical).

### R2. `iface_gen_method_gen` — generated element mis-classified as input

- **Corpus**: validation-suite 233.
- **Was**: the reconstruction's drop-list turned the **generated** interface
  element 3 into an input element and dropped the genuine input
  `element_geometry` records, so the target `element_geometry 3 = 10` did not
  exist ("Error in calculation").
- **Resolution**: remove the input `element 3` line and restore the input
  `element_geometry 1 20` / `2 30`; the generator now produces element 3 with
  `element_geometry 3 10`, so the target reads the **generated** record (the run's
  `element_geometry` set is byte-identical to the archive).
- **Commit**: `fbcf144`.
- **Evidence**: `CORPUS-STALE-TARGETS.md` §2.

### R3. `elmethod` / `numit` — lost `-follow_material` reference geometry + initial truss force

- **Corpus**: validation-suite 233.
- **Was**: the reconstruction used the dump `node` (the **deformed** geometry
  under `-follow_material`) as input, and the drop-list removed the initial
  `element_truss_force` self-stress (the dump keeps only the relaxed final value
  ≈0). The secondary velocity targets were therefore unrecoverable.
- **Resolution**: restore the reference geometry (`node_start_refined` is the
  input reference; `top.cc:269`) and derive the initial forces from the
  documented relaxation law `F_init = k·(L_init − L_final)`, `k = E·A/L_init`
  (elmethod 0.123 / 0.246, numit 0.123). The reconstructed inputs reproduce the
  archived `.dbs` byte-for-byte (state and prints).
- **Commit**: `c8ee32e`.
- **Evidence**: `CORPUS-STALE-TARGETS.md` §"`-follow_material` reference geometry
  + initial truss force" (validation-suite 234: 231/3 → 233/1; the 233 real
  tests 233/233; sfnet unchanged).

### R4. Post-solve stress pass cost (old C4) — hoisted out of the loop

- **Corpus**: all models with `materi_stress` + the absolute recovery (e.g. the
  delete/mpc family).
- **Was**: the post-solve stress sync (`dae06a3`) ran a stress-only
  `element_loop()` after every equilibrium iteration, roughly doubling the
  element-assembly cost for those models.
- **Resolution**: the pass was hoisted out of the equilibrium loop (`3957809`,
  docs `e0a3997`): it now runs once per step and assigns the stress dofs directly
  via a stress-only variant of `parallel_new_dof_diagonal`. Measured speedups:
  `large3` 1.23–1.25×, `dynamic8` 1.12–1.21×, `delete3` ~1.14× (small tests are
  startup-dominated). This is improvement item (a) in `IMPROVEMENT-ANALYSIS.md`,
  now marked **DONE**.
- **Acceptance**: mpc5 `sigyy=0.5` at k=1,2,3,5, delete2 `sigxx=0.8`, delete3
  `sigyy=0.6` unchanged; sfnet 197 and validation-suite 233 classification diff
  0; internal suite 16/16. Two non-target derived-dof changes are documented
  (`contact_heatgen` `xsigyy`; `mpc8` `sigxx`, which now matches the Professional
  better).
- **Evidence**: `manual-developer/nodal_stress_post_solve_sync.md`
  §Cost / §Verification.

### R5. Damage stress compounding — legacy 2014 defect (2026-09-29)

- **Corpus**: sfnet `damage1` (target checks `dam` only, so it never failed);
  found by the keyword convergence cycle (`materi_damage_apply` measurements
  of 2026-09-27: `dam` identical to the Professional, `sigxx` not).
- **Was**: `set_stress()` multiplied the carried, already damaged stress by
  (1-d) every step (`sigma_n = (1-d_n)(sigma_{n-1} + C deps)`): the `damage1`
  bar gave `sigxx` 1436.83 instead of (1-d)E eps = 85598.32 (Professional
  85598.32); at constant strain the stress decayed to ~0. Present in
  `Sources-2014.zip`.
- **Resolution**: the effective (undamaged) stress is carried per
  integration point in the internal record
  `element_intpnt_materi_stress_effective`; the damaged stress is
  (1-d) x effective, never divided. Yield test on the damaged stress and
  gated-step semantics measured on the Professional 02-08-2026. Commit
  pending.
- **Acceptance**: `validation-suite/test-2026/damage_stress1.dat`,
  `damage_stress_hold1.dat`, `damage_stress_3d1.dat` pass on the GNU and the
  Professional; validation-suite 233 and sfnet 187/10 classification diff 0;
  the only real `.dbs` change is sfnet `damage1` `sigxx`.
- **Evidence**: `manual-developer/materi_damage_effective_stress.md`.

### R6. `-updated_linear` / `-updated_area` stress rotation and compressibility (2026-09-29)

- **Corpus**: audit cluster K01 (13 tests: Pro `interac1`, `matrix7`; sfnet
  `druckpr1`, `interac1`, `matrix7`, `mohrcou1`, `tria6`; validation-suite
  `mdp_shear`, `mesh_act_grav`, `mshf`, `mshf_nof`, `mstrain_druckprag`,
  `mstrain_druckprag_elast`, all on the converted `-updated_linear` input).
- **Was**: (1) `materi.cc` rotated the `-updated_linear`/`-updated_area`
  stress by `inc_rot`, which `set_deften_etc()` overwrites with the real
  incremental rotation (added by `22801c3`/`a1d9d12`); simple shear
  `druckpr1` gave sigxy 0.0693, sigxx 0.0924 instead of 0.11547 / 0.
  (2) `stress.cc` compressibility of the linear memories used
  `tr(new_epe-old_epe)` with `old_epe` = 0 unless `materi_strain_elasti` is
  initialized, so the total strain was added every step (legacy 2014):
  `interac1` -0.01 instead of -0.005.
- **Resolution**: the linear memories are removed from the rotation list
  (manual Pro 2.2.1: no rotation step); the volumetric increment is
  `tr(inc_epe)/compressibility`. Commit pending.
- **Acceptance**: `validation-suite/test-2026/updated_linear_compress1.dat`,
  `updated_linear_shear1.dat` pass on the GNU and the Professional (rc=1 on the
  pre-fix binary); 9 of the 13 K01 tests now match the Professional in every
  variable, the 4 `-fixed_in_space` validation-suite tests keep a residual
  identical to `-updated_without_rotation`; validation-suite 233 and sfnet
  187/10 classification diff 0, no physical `.dbs` change.
- **Evidence**: `manual-developer/group_materi_memory_updated_linear.md` §5,
  `AUDIT-ALL-VARIABLES-2026.md` K01.

### R7. Static phreatic level without a pressure unknown (2026-09-29)

- **Corpus**: audit cluster K22, Pro `ground18` (`ground8` of the same
  cluster has a pressure unknown and is not affected; see the audit).
- **Was**: with `groundflow_phreatic_level_static -yes` and no
  `groundflow_pressure`, `get_materi_density()` chose the wet density for
  every element (the `NODE_PHREATICLEVEL` record needs the pressure
  unknown) and `materi()` added no pore pressure to the skeleton. Bottom
  sigyy -1995 instead of -1247.5 (analytic -1250).
- **Resolution**: for groups with `-groundflow`
  (`groundflow_static_pressure_group()`), the static total pressure of
  `groundflow_phreatic_coord()` drives the wet/dry choice and acts on the
  skeleton (scaled by `groundflow_pressure_factor`). Commit pending.
- **Acceptance**: `ground18` equals the Professional in every variable;
  `validation-suite/test-2026/phreatic_static1.dat` and
  `phreatic_static_submerged1.dat` pass on the GNU and the Professional
  (rc=1 on the pre-fix binary); validation-suite, sfnet and test-2026
  `.dbs` identical to the pre-fix binary apart from the 2 new tests.
- **Evidence**: `manual-developer/groundflow_phreatic_level_static.md`,
  `AUDIT-ALL-VARIABLES-2026.md` §5.2. Side effect: C10.

### R8. Control records honored for index 0 only (2026-09-29)

- **Corpus**: audit cluster K04, Pro `elasti7` (`control_materi_elasti_k0
  20 -yes`).
- **Was**: `db_active_index( CONTROL_xxx, 0, ... )` tests that exact
  index. `control_materi_elasti_k0` was ignored for any other control
  index (elasti7 sigxx 0 instead of -0.05);
  `control_materi_plasti_hardsoil_gammap_initial` with a non-zero index
  aborted ("Data is allocated in a parallel loop");
  `control_solver_bicg_stop` followed neither the manual nor the
  Professional for mixed indices / with `solver_bicg_stop`.
- **Resolution**: K0 switch read at ICONTROL, and (measured on the
  Professional) applied without a poisson record and to hardsoil
  nu50/nuur; gammap pre-allocation when the record exists at any index;
  bicg stop = control record index 0 overruling the global one (the
  Professional's measured semantics, `47b6fda`), then extended to the
  combined precedence of C13 (commit pending).
- **Acceptance**: elasti7 equals the Professional in every variable;
  `validation-suite/test-2026/elasti_k0_index1`, `elasti_k0_nopoisson1`,
  `elasti_k0_hardsoil1`, `hardsoil_gammap_initial_index1` pass (rc=1 on
  the pre-fix binary); 9/9 bicg-stop probes matched the Professional
  with `47b6fda` (5/9 before; 6/9 with the C13 combined precedence);
  validation-suite, sfnet and test-2026 `.dbs` identical to the pre-fix
  binary apart from the 4 new tests.
- **Evidence**: `manual-developer/group_materi_elasti_k0.md`,
  `solve_iterative_bicg.md`, `AUDIT-ALL-VARIABLES-2026.md` §5.3. Found on
  the way: C11, C12.

### R9. `_direct_normal_automatic` inert in 2D — `element_normal` (2026-09-30)

- **Corpus**: audit cluster K04, Pro `other/tension_direct1`
  (`group_materi_plasti_tension_direct_normal_automatic`).
- **Was**: `materi()` took the "automatic" plane normal as the cross
  product of two element edges: (0,0,+-1) in 2D, a direction with no
  stress, so the cut-off never acted (`tension_direct1`: GNU elastic,
  `disy` 0.42 at the middle nodes, `sigyy` 0.667; Pro `disy` 0, all
  stresses 0). The loop also read `NODE` coordinate 2 of 2D nodes (out of
  bounds). The GNU had no `element_normal` record (keyword gap).
- **Resolution**: new generated record `element_normal` (manual 6.445)
  with the rule measured on the Professional (boundary side = side not
  shared with another isoparametric element; node normal = normalized sum
  of the consistent side normals; element normal = normalized sum of its
  node normals; regenerated after mesh changes only), computed in
  `iteration_start()` before the parallel element loop and read by
  `materi()` (commit pending).
- **Acceptance**: `tension_direct1` equals the Professional in every
  variable except `eptyy` at the middle nodes (0.667 vs 0.5, nodal strain
  recovery, audit K11); element_normal equal to the Professional to print
  precision on ~30 probe meshes; 5 new `validation-suite/test-2026/
  element_normal*` inputs pass on the Professional and the GNU (the
  pre-fix binary rejects `-element_normal`; without those targets it
  fails `disy`); validation-suite 233 PASS, sfnet 187 PASS / 10 RUNFAIL
  and test-2026 23/23 unchanged (28/28 with the new inputs), all
  `.dbs` identical to the pre-fix binary except `materi_direct_auto`
  (single hex8: new zero `element_normal` record, results changed at
  round-off level only, 2.8e-14).
- **Evidence**: `manual-developer/element_normal.md`. Residuals: C14.

### R10. Default material memory: `-updated` instead of the Professional `-updated_linear` (old C16, 2026-10-04)

- **Was**: with no `group_materi_memory` record the GNU used `-updated`
  (hard-coded at five places), the Professional behaves as
  `-updated_linear`: single quad4 simple shear gamma 1, GNU
  0.2311 / 0.4207 / -0.2286 vs Pro 0 / 0.5 / 0 (sxx / sxy / syy); the
  2-quad4 column of C14 sheared by 1 was nonlinear in the GNU only.
- **Resolution**: `materi_memory_default()` (general.cc) gives
  `-updated_linear`, or `-updated` with `control_materi_updated_apply
  -yes` at the current timestep (manual 6.152); assigned before every
  memory read (materi, set_deften_etc, stress, elem, tendon); the
  reference-geometry gradient stays limited to an explicit
  `-updated_linear` (`d8117aa`).
- **Acceptance**: `validation-suite/test-2026/materi_memory_default_shear1`
  (= Pro; rc=1 on the pre-fix binary), `materi_memory_default_stretch1`
  (= Pro 0.695653), `materi_memory_default_updated_apply1` (manual
  semantics, see C21); validation-suite 233/1, test-2026, sfnet 187
  unchanged; Pro `other` mohr_coul_direct6/8 now pass, `hypo12` becomes
  a marginal miss (C22).
- **Evidence**: `manual-developer/group_materi_memory_default.md`,
  `control_materi_updated_apply.md`. Residuals: C18-C22.

---

## B. sfnet 197 — the 9 RUNFAIL (all out-of-scope)

### B1. `examp6` — contact + adaptive remesh node-numbering regression

- **Corpus**: sfnet 197.
- **Symptom**: target `node 259 x = 0.51 ±0.02`; GNU-current `0.609375`
  (299 nodes). The reference GNU 2014 and the repo's first modern build give
  `0.5078125` (274 nodes, rc=0).
- **Root cause**: regression from the `ELEMENT_DOF` stress/history slot fix
  `8f16dc0` (2026-08-29), compounded by the `inertia_apply` default
  `0560dcf`. `8f16dc0` changed the remembered integration-point stress ranges
  and the constitutive write index for `nder>1`; examp6 uses `derivatives`
  (2D ⇒ `nder=4`), so the mapped material state across the local
  refinement/restart changed, which changed the equation residue and therefore
  which elements the residue-driven `control_mesh_refine_locally` selected.
  The target coordinate still exists in the mesh, but at node 274 (default) /
  262 (`inertia_apply -yes`).
- **Evidence**: `DIAG-CONTACT-REFINEMENT.md` §1; bisect first-bad `8f16dc0`.
- **Severity**: **MEDIUM** — the target is tied to a solver-sensitive
  adaptive-remesh numbering; the physics path is a correct fix and reverting it
  would regress the `gforce7/10/13` and `qsri3d` family.
- **Proposed fix**: recalibrate the target to the corrected remesh numbering,
  or a dedicated adaptive-remesh work unit that makes the node numbering
  residue-stable. No minimal source fix exists.
- **Effort**: **cheap** (target) / **large** (adaptive-remesh work unit).

### B2. `examp14` — stale target + honest-solver runtime

- **Corpus**: sfnet 197.
- **Symptom**: 45 s timeout; with a 120 s budget the value is `9.71887788387`
  (current) / `9.61289226266` (2014), while the stored target is `15.131 ±1`
  (no GNU build reaches it, difference >50%). Runtime 47.6 s vs the 45 s slot.
- **Root cause**: the stored target is the original 2000-08-21 value and
  matches no available binary (the Professional tree has no `examp14.dat` and
  rejects the GNU `options_mesh` dialect). The extra runtime is the honest
  Bi-CG converging to machine zero (80 iterations, error 9.9e-17) instead of
  the removed 2014 stagnation exit (8 iterations, error 9.4e-13).
- **Evidence**: `CORPUS-SOLVER-TIMEOUT.md` §2, §6.1.
- **Severity**: **LOW** — a reference/target mismatch, not a fixable source
  defect.
- **Proposed fix**: recalibrate the target to the honest GNU value and extend
  the corpus timeout for this test (harness budget), or drop it from the
  regression set.
- **Effort**: **cheap**.

### B3. `examp17` — stale target (contact + p-refinement)

- **Corpus**: sfnet 197.
- **Symptom**: target `post_point_dof 0 -sigyy = -20506 ±100`; GNU-current
  `-21103.3`, GNU 2014 `-21330.5`.
- **Root cause**: the stored target is met by no available binary. The
  Professional cannot run the corpus file (GNU dialect); its own equivalent
  `test/validation/validation_10.dat` uses a different model and a 100×
  looser target (`-20538.3 ±1e4`), giving `-21031.0` (no p) / `-20699.7`
  (with p).
- **Evidence**: `DIAG-CONTACT-REFINEMENT.md` §2.
- **Severity**: **LOW** — stale target, not a source defect.
- **Proposed fix**: align the input/target with the Professional
  `validation_10` semantics (a separate convergence work unit), or recalibrate
  to the honest GNU value.
- **Effort**: **medium**.

### B4. `examp24` — O(N²) mesh mapping makes the test intrinsically long

- **Corpus**: sfnet 197.
- **Symptom**: 45 s timeout; 150 `control_repeat` remesh+step cycles. One
  cycle is 16.09 s (current) / 7.79 s (2014); 150 cycles ≈ 40 min / 20 min.
- **Root cause**: the dominant cost is the mesh-to-mesh **mapping**
  (`new_mesh_version` → `parallel_map_element/node` → `point_el` →
  `tetrahedron_volume` → `matrix_ab`), a brute-force O(N_new × N_old)
  point-in-element scan. Instrumented: 29,065,922 `point_el` calls per cycle,
  0.08% found; `matrix_ab` is 51.4% of sampled runtime. `map.cc` is
  byte-identical to the 2014 reference, so the O(N²) is inherited, not a
  regression.
- **Evidence**: `CORPUS-SOLVER-TIMEOUT.md` §2, §6.2.
- **Severity**: **LOW** — a performance limitation with no target; no minimal
  source fix.
- **Proposed fix**: a spatial index / bounding-box pre-filter in the mapping
  path. Even then the test remains ~20 min, so the corpus slot must also
  change.
- **Effort**: **large**.

### B5. `ho_mech2` — stale p-refinement instability input/target

- **Corpus**: sfnet 197.
- **Symptom**: target `vely = -641`; GNU-current `-1238`.
- **Root cause**: stale input/target. The Professional's own copy
  (`test/other/tensile_instability.dat`) uses `-updated_linear` + a different
  timestep control and a different target (`-1975`). The 2014 GNU aborts
  (`number_of_intagration_points at least 9`); the direct solver also gives
  `-1238`. No binary reaches `-641`.
- **Evidence**: `CORPUS-PARSE-INIT-EXAMP2.md`.
- **Severity**: **LOW**.
- **Proposed fix**: align the input/target with the Professional's
  `tensile_instability.dat` variant.
- **Effort**: **medium**.

### B6. `ho_mech3` — model divergence present in the reference GNU

- **Corpus**: sfnet 197.
- **Symptom**: 45 s timeout; the run diverges (`post_node_result 1 -vely` ≈
  1e15 vs target 22000).
- **Root cause**: both GNU builds diverge in the same `so.cc` guard
  (`EPS_SOL 1e15`) inside `hyperelasticity`/`plasti_rule`. The 40-cycle model
  is unstable in both; the target is transiently reached at 2-5 cycles. The
  Professional cannot run the model (rejects `bounda_unknown`; no
  `control_mesh_new_mesh` test in its tree).
- **Evidence**: `CORPUS-SOLVER-TIMEOUT.md` §2, §6.3.
- **Severity**: **LOW** — model/numerics limit, not a GNU regression.
- **Proposed fix**: none minimal; a stabilisation/model fix would be a
  dedicated work unit.
- **Effort**: **large**.

### B7. `ho_mech4` — missing `materi_strain_plasti`; stale input/target

- **Corpus**: sfnet 197.
- **Symptom**: target `kap = 0.0747218`; GNU-current `kap=0` (init abort);
  with the missing record added, `kap=0.06274`.
- **Root cause**: the initia declares `materi_strain_total` but not the
  `materi_strain_plasti` required by `group_materi_plasti_vonmises` (manual
  Pro 6.750). The Professional's own copy (`test/other/cylinder_buckle.dat`)
  declares it and uses a different mesh + target (`node_dof 1010 -velx =
  0.0118`). The sfnet target is not the Professional's target.
- **Evidence**: `CORPUS-PARSE-INIT-EXAMP2.md`.
- **Severity**: **LOW** — stale input/target.
- **Proposed fix**: align the input/target with the Professional's
  `cylinder_buckle.dat`.
- **Effort**: **medium**.

### B8. `petsc1` — PETSc not a build option in this environment

- **Corpus**: sfnet 197.
- **Symptom**: `Error: set PETSC_USE to 1 in tnpetsc.h`.
- **Root cause**: `tnpetsc.h` defines `PETSC_USE 0`; `so_petsc.c` is entirely
  inside `#if PETSC_USE`; no PETSc headers/libraries are installed and the
  Professional 02-08-2026 tree ships no `tnpetsc.h` either.
- **Evidence**: `CORPUS-SOLVER-TIMEOUT.md` §3, §6.4.
- **Severity**: **LOW** — build infrastructure, not a source defect.
- **Proposed fix**: install MPI+PETSc and wire the build (note: enabling PETSc
  forces SuperLU off, `so.cc:220-221`, losing the direct solver used by the
  interface/patch tests). Cost >> benefit for one 3-node 1-D heat test.
- **Effort**: **large**.

### B9. `taylor3` — stale input + model/solver regression

- **Corpus**: sfnet 197.
- **Symptom**: target `kap = 2.58442`; GNU-current stalls at the minimal
  timestep (`kap=0`); GNU 2014 gives `2.56391` (fails).
- **Root cause**: stale input/target. The Professional's own copy uses
  `inertia_apply`, `mesh -follow_material`, `group_materi_memory -updated`, a
  3-value timestep control and `_stop -no`; the Professional rejects the GNU
  input (`options_stabilization`). No binary reaches the target.
- **Evidence**: `CORPUS-PARSE-INIT-EXAMP2.md`.
- **Severity**: **LOW**.
- **Proposed fix**: align the input with the Professional's variant, then
  re-evaluate the residual as a solver/model work unit.
- **Effort**: **medium**.

### B10. sfnet targets calibrated on the legacy `convection_apply` default -yes (K03, 2026-10-06)

- **Corpus**: sfnet 197 (external, not edited). Since K03 the GNU default
  `convection_apply` is `-no` (manual Professional 6.395 and the
  Professional 02-08-2026). Five inputs without the record relied on the
  legacy GNU default `-yes` and now miss their targets (sfnet 187/10 ->
  182/15):

  | test | target | GNU old (-yes) | GNU new (-no) | Professional |
  |---|---|---|---|---|
  | `condif2` | node 2 temp 1 +/- 0.02 | 1 | 0.5 | 0.5 |
  | `condif3` | node 2 temp 1 +/- 0.02 | 1 | 0.5 | 0.5 |
  | `examp1` | post-point velx -0.041 +/- 0.1 | -0.04873 | 0.1391 | 0.1391 |
  | `incnav5` | post-point velx 0.389082 +/- 0.01 | 0.3891 | 0.4919 | 0.4919 |
  | `strconv1` | node 5 sigxx 1 +/- 0.01 | 1 | 0 | RUNFAIL (pardiso error) |

- **Root cause**: stale targets (input dialect), not a GNU defect: the
  new values equal the Professional on the same (converted) input.
- **Proposed fix**: none in the GNU; a versioned copy of these inputs
  would need an explicit `convection_apply -yes` to keep testing the
  convection.
- **Severity**: **LOW**. **Effort**: **cheap**.

### B11. sfnet `incnav3` target predates the `expansion_volume` inertia density (K04, 2026-10-06)

- **Corpus**: sfnet 197 (external, not edited). Since F1 T3 the GNU uses
  rho = rho0*(1 - beta*T) of `group_materi_expansion_volume` in the
  momentum inertia too (Professional 02-08-2026; its own `incnav3` states
  "the density gets value rho = initial_rho*{1-expansion*T} = 0.5 ... the
  total rate is -5"). The legacy sfnet copy keeps the old comment and
  target `tvely` -2.5 +/- 0.01: GNU old -2.4976, GNU new -4.9904, the
  Professional on the converted sfnet input -4.9904 (sfnet 181/16 ->
  180/17 with the audit runner).
- **Root cause**: stale target, not a GNU defect.
- **Proposed fix**: none in the GNU.
- **Severity**: **LOW**. **Effort**: **cheap**.

### B13. sfnet targets calibrated on the old GNU mesh/memory defaults (F3, 2026-10-07)

- **Corpus**: sfnet 197 (external, not edited). Since F3 the default mesh
  is `-fixed_in_space` (manual 6.806) and the default memory with
  `materi_displacement` is `-total_linear` (manual 6.797), as the
  Professional. Eight inputs without a `mesh` record relied on the legacy
  GNU moving mesh (`-follow_material`) and now miss their targets (sfnet
  exit-code pass 180/181 -> 173):

  | test | target | GNU old | GNU new | Professional 02-08-2026 |
  |---|---|---|---|---|
  | `axisym3` (`-updated`) | node 2 velx -1e-3 +/- 1e-4 | -0.000999951 | -0.00999896 | does not run the sfnet input (contact target record) |
  | `axisym4` (`-updated`) | node 2 velx -1e-3 | -0.000999986 | -0.00199997 | does not run the sfnet input |
  | `contac1` | node 2 x 1.01 +/- 1e-3 | 1.01 | 0.99 | does not run the sfnet input (contact target record) |
  | `contac2` | node 2 x 1 | 1 | 0.99 | 0.99 (= GNU new) |
  | `examp12` (`-updated`) | post-point sigxy -1.13905e-4 | -1.18897e-4 | -6.80231e-4 | does not run the sfnet input |
  | `temp4` (`-empty` group, condif) | node 1 xtemp 0.5 | 0.5 | 1 | 0 (divergence, see DIVERGENCES-PROFESSIONAL.md) |
  | `temp5` (`materi_displacement`, condif) | node 1 xtemp 0.5 | 0.5 | 1 | 0 (same) |
  | `visc_pl2` (`-updated_without_rotation`) | node 1 sigxx 1.66666 | 1.66666 | 2 | does not run the sfnet input (parse error) |

- **Root cause**: stale targets of the external corpus (they encode the
  legacy GNU `-follow_material` default); not a GNU defect, except the
  `temp4/5` condif difference listed in DIVERGENCES-PROFESSIONAL.md.
- **Proposed fix**: none in the GNU; a versioned copy would need an
  explicit `mesh -follow_material ...` record.
- **Severity**: **LOW**. **Effort**: **cheap**.

---

## C. Additional residuals

### B12. validation-suite `node_force_inertia` target calibrated on the reversed `node_force` sign (C28, 2026-10-06) — RESOLVED 2026-10-06

- **Corpus**: validation-suite test-2014 (run by `build_safe.sh`). Its
  `node_inertia 3` target -19.83 +/- 0.2 was calibrated when the GNU
  applied `node_force` in the opposite direction (the -30 force pushed
  node 3 UP: -19.75). Since C28 the force pushes down, adding to the
  gravity of the node mass and of the element: -20.109 (outside the
  tolerance by 0.08). The Professional cannot run the input (`node_mass`
  rejected, `node_inertia` wants one value per dof).
- **Root cause**: stale target, not a GNU defect.
- **Proposed fix**: target -20.11 +/- 0.2 (one line in
  `validation-suite/test-2014/node_force_inertia.dat`).
- **Resolution** (project owner go-ahead, 2026-10-06): target retargeted
  to -20.11 +/- 0.2 with a comment in the input, together with the C28
  fix on `dev/f1-cheap-bugs` (F1 T6).
- **Severity**: **LOW**. **Effort**: **cheap**.

### B14. sfnet targets calibrated on the old GNU `-size_tot`/`-size_dev` (F4 T6a, 2026-10-11)

- **Corpus**: sfnet 197 (external, not edited). Since F4 T6a `-size_tot`
  and `-size_dev` are sqrt(A:A) / sqrt(B:B) as the manual (post_calcul)
  and the Professional; the GNU used sqrt(0.5 ...). Four sfnet inputs
  have targets calibrated on the old value and now exit 1:
  `calcul1`, `calcul2` (sizes 0.707 / 0.577; the Professional's own
  `calcul1` targets 1 / 0.816497), `ctvmis2` (`-sizedev` 0.11547; the
  Professional's own `ctvmis2` targets 0.163299 = 0.11547 sqrt(2), which
  the GNU now meets: Pro `other` ctvmis2 exit 1 -> 0) and `examp18`
  (`post_point_dof_calcul` 6.07489e7; GNU node values 1.514e8 ->
  2.142e8, ratio sqrt(2)). Other size values change by sqrt(2) with no
  exit-code change (`examp3`, `examp10`, `examp14`, `examp17`,
  `examp23`).
- **Proposed fix**: retarget the four inputs (x sqrt(2)) when the sfnet
  corpus is next edited.
- **Severity**: LOW. **Effort**: cheap.

### B15. sfnet `truss2`/`truss3` targets calibrated on the legacy GNU truss memory default (F4 T6d, 2026-10-11)

- **Corpus**: sfnet 197 (external, not edited). Both inputs have no
  `group_truss_memory` record and target a reaction (0.1, 0.1) that only a
  geometrically nonlinear truss gives. Since F4 T6d the default is the
  Professional one (`-updated_linear` without, `-total_linear` with
  `materi_displacement`): GNU reaction (0.1, 0), exit 0 -> 1. The
  Professional 02-08-2026 runs the same inputs (dialect-translated) to
  (0.1, 0) too; its own `other/truss2`/`truss3` add `group_truss_memory
  0 -updated` (and `truss3` replaces `materi_displacement` by
  `materi_velocity_integrated`).
- **Proposed fix**: adapt both inputs as the Professional copies when the
  sfnet corpus is next edited.
- **Severity**: LOW. **Effort**: cheap.

### C1. `pruebavel` — dead/malformed test (the 1 PARSE)

- **Corpus**: sfnet 197.
- **Symptom**: unknown-record error (`group_condif_flow 0 -velx -vely` uses
  dof names; manual Pro 6.600 is numeric `beta1 beta2 beta3`). The
  validation-suite copy is 0 bytes.
- **Root cause**: a dead/malformed test. GNU, GNU 2014 and the Professional all
  reject the record; no dof-name variant is documented.
- **Evidence**: `CORPUS-PARSE-INIT-EXAMP2.md`.
- **Severity**: **LOW** — test artifact.
- **Proposed fix**: none; keep it classified as `PARSE` (it is not a
  regression). Optionally remove it from the regression set.
- **Effort**: **cheap**.

### C2. `groundflow_interface` — `gvelx` sign convention residual

- **Corpus**: validation-suite 233 (the test itself now PASSes on pressure).
- **Symptom**: the `gvelx` dof sign differs from the archived Professional dump
  (`-1.9996159999e-06` vs `+1.99960004065e-06`); the pressure field is
  identical and matches the archive to 5.4e-8.
- **Root cause**: not established. The pressure target is met; the residual is
  a velocity-sign convention issue unrelated to the pressure and to the
  reconstruction. (The main groundflow velocity sign fix `0838457` is
  separate and closed.)
- **Evidence**: `CORPUS-RECONSTRUCTION.md` fifth pass.
- **Severity**: **LOW** — no target depends on it.
- **Proposed fix**: investigate the sign convention of `groundflow_velocity`
  in the interface-coupled path against the Professional dump.
- **Effort**: **medium**.

### C3. `gforce10` — 3D section face/arm sign limitation

- **Corpus**: sfnet 197 (target recalibrated to the honest GNU value).
- **Symptom**: `node_dof_calcul 5` value 9 (`mom1y_sig`): GNU returns `+5000`
  where the Professional `force10.dat` returns `-5000`; the magnitudes and
  `nory = -12.34` agree.
- **Root cause**: the GNU 3D section post-processing has a face/arm sign
  issue. Documented in `DIAG-SOLVE-MIXTO.md` §13.4 ("a remaining 3D face/arm
  issue of the section post-processing, NOT the scheme; the N is exact"). The
  target was recalibrated to `+5000` so the test guards the honest GNU value.
- **Evidence**: `CORPUS-STALE-TARGETS.md` §3.
- **Severity**: **MEDIUM** — a post-processing correctness divergence a user
  reading section moments would hit.
- **Proposed fix**: correct the 3D section face/arm sign in the section
  post-processing, then restore the Professional-aligned target (`-5000`).
- **Effort**: **medium** (post-processing subsystem).
- **Status (F4 T6a, 2026-10-11)**: the sign is the
  `post_calcul_materi_stress_force_plot_switch` convention and the frame
  needs `post_calcul_materi_stress_force_thickness_switch -yes`; the
  test-2014 `gforce10.dat` had lost both records of the Professional
  `force10.dat`. Restored, target back to `-5000`: GNU node 5 nory
  -12.34, shey -100, mom1y -5000 = Professional exactly. The
  thickness_switch default is `-no` in the manual (6.915: "normally ...
  the shortest element direction"), in the Professional (force10 with
  `-no` = without the record) and in the GNU: no change. Remaining
  difference, OWNER DECISION pending: a SQUARE end face without `-yes`
  (no shortest direction; the manual is silent) - the Professional takes
  a degenerate frame (vector components 0, the moment in the mom2 slot:
  force10 without the switch gives nors -12.34, mom2s 5000, all vectors
  0), the GNU breaks the tie with the reference point (nory -12.34, mom1y
  5000). Owner decision 2026-10-11: the GNU keeps its tie-break
  (deliberate divergence, class A in DIVERGENCES-PROFESSIONAL.md). Closed.
  The sfnet copy of `gforce10` (on omen, not versioned) keeps the
  old input and target.

### C5. `contact_target_geometry` alias only reads index 0

- **Corpus**: not exercised by the corpus (all corpus tests use the GNU name
  `contact_geometry`); it affects the Professional-dialect `validation_10`.
- **Symptom**: the GNU silently drops a `contact_target_geometry 1 ...` record
  (index 1) and returns ~0 for `validation_10`.
- **Root cause**: the alias added in `contact.cc` copies only index 0
  (`db_active_index(CONTACT_TARGET_GEOMETRY, 0, ...)`).
- **Evidence**: `DIAG-CONTACT-REFINEMENT.md` "Side finding".
- **Severity**: **MEDIUM** — silently drops a boundary condition on
  Professional-dialect input.
- **Proposed fix**: move the alias out of the parallel region
  (`parallel_contact` runs through `parallel_sys_routine`; the `PUT` allocates)
  and handle all indices.
- **Effort**: **medium**.
- **Status**: **DONE 2026-10-07 (F2 T3a)**. Broader than first documented
  (audit K05): index 0 aborted with "Data is allocated in a parallel loop"
  (the alias PUT ran inside `parallel_contact`), index > 0 was dropped. Fix:
  read-only `contact_geometry_get()` in `contact.cc` (both names, every
  index, GNU name wins on the same index); no record is created. Tests
  `contact_target_geometry_index1/_index0/_switch1`; validation_10 sigyy
  -21232 (Pro -21031). Element-target contact (`contact_target_element_group`)
  had a second, separate root cause (no mesh boundary): fixed in F2 T3b,
  residuals in C40.

### C6. `iface_gen*` — reconstructed inputs do not exercise interface generation — RESOLVED 2026-10-11 (F4 T6a)

- **Corpus**: validation-suite 233 (`iface_gen`, `iface_gen_method`,
  `iface_gen_geom`, `iface_gen_quad9`).
- **Symptom**: the same drop-list defect that produced R2 turned the
  **generated** interface element(s) into input elements in these four tests as
  well (`iface_gen` element 3, `iface_gen_method` element 3, `iface_gen_geom`
  element 3, `iface_gen_quad9` elements 3 and 4). Their targets are `node_dof`
  values, so they pass either way, but the reconstructed input does not exercise
  the generation feature.
- **Root cause**: `element_geometry` was classified as a run-produced record and
  dropped; the generator's output element was kept as an input element.
- **Evidence**: `CORPUS-STALE-TARGETS.md` §2 "Family note".
- **Severity**: **MEDIUM** — a corpus-coverage gap: a future regression in
  interface generation would not be caught by these tests.
- **Proposed fix**: apply the R2 correction to the four siblings (remove the
  generated input element, restore the input `element_geometry` tags) and confirm
  the targets still pass; five currently-passing tests are touched, so it is a
  separate work unit with a full A/B.
- **Effort**: **cheap**.
- **Status**: **RESOLVED 2026-10-11 (F4 T6a)**: the generated elements
  (`iface_gen`, `iface_gen_method`, `iface_gen_geom`: element 3;
  `iface_gen_quad9`: elements 3 and 4) and their `element_group ... 10`
  records were removed from the inputs; `iface_gen_method` (selection
  `-element_geometry`) got back its input tags `element_geometry 1 20 /
  2 30`. Each run generates the same elements (`element 3 -quad4 2 4 5
  7`, quad9: `3 -quad4 6 9 13 16`, `4 -quad4 3 6 10 13`, `element_group
  3/4 10`), rc 0, targets met. build_safe checks that the input has no
  element 3 and the run produces it.

### C7. Mazars damage: gamma with shear, and damage + early plasticity

- **Corpus**: none (own hex8 / bar probes, `materi_damage_effective_stress.md`).
- **Symptom**: (a) with a compressive `sigxx` and a same-sign `sigxy`/`sigxz`
  the Mazars `dam` differs from the Professional 02-08-2026 (e.g. 0.10809 vs
  0.18752; with the opposite shear sign 0.26351 vs 0.39352); (b) a bar with
  `group_materi_plasti_vonmises` yielding before the damage threshold gives a
  different plastic strain path (`eppxx` 1.879e-4 vs 1.368e-4).
- **Root cause**: (a) `damage_mazars()` builds gamma from xx, xy, xz (not
  objective); the Professional is not objective either and matches none of
  the candidate forms tried (GNU form, principal-stress gamma, gamma = 1).
  (b) not established (iteration order of plasticity and damage); both codes
  satisfy sigma = (1-d) C eps_e for their own strains.
- **Severity**: **MEDIUM** — damage in shear-compression states and
  damage+plasticity models diverge from the Professional.
- **Proposed fix**: identify the Professional gamma with more probes (vary
  each stress component separately), then (b) with a trace of the plastic
  iterations.
- **Effort**: **medium**.

### C8. Default `mesh`: GNU `-follow_material`, Professional `-fixed_in_space` — RESOLVED 2026-10-07 (F3 T2)

- **Resolution**: top.cc default `-FIXED_IN_SPACE` (manual 6.806); tests
  `mesh_default_fixed_in_space1/2` (= Pro). Corpus consequences in B13.

- **Corpus**: none directly; explains the `interac1` `-updated` residual
  (GNU -0.0050117, Pro -0.005) and why the `c945833` reference gradient of
  `-updated_linear` matches the Professional only without a `mesh` record.
- **Symptom**: with `materi_velocity` and no `materi_displacement` the GNU
  moves the nodes with the material by default (`top.cc`, `options_mesh`
  default `-FOLLOW_MATERIAL`); the Professional manual 6.806 says "Default
  each specifier is set to -fixed_in_space", and the Professional keeps the
  nodes in place. 1D bar stretched 50% (dt 0.01): `-updated` Pro 0.5 / GNU
  0.40714 without a mesh record, both 0.40714 with `mesh -follow_material`.
  With `mesh -follow_material` the Professional `-updated_linear` gives the
  current-geometry value 0.40714 while the GNU (reference gradient, `c945833`)
  gives 0.5.
- **Root cause**: different mesh default; the GNU `-updated_linear` reference
  gradient emulates the Professional's fixed mesh for that one memory only.
- **Severity**: **MEDIUM** — any `-updated`/velocity-only model without a
  `mesh` record differs at large strain; an explicit `mesh -follow_material`
  with `-updated_linear` differs too.
- **Proposed fix**: align the default mesh with the Professional (then drop
  the reference-gradient special case of `-updated_linear`); needs a full
  corpus A/B, because many GNU tests rely on the moving mesh.
- **Effort**: **large**.

### C9. `-total_linear` + compressibility: the Professional is time-step dependent

- **Corpus**: none (no corpus test combines them); probe = Pro `interac1` with
  `-total_linear` and `materi_displacement`.
- **Symptom**: analytic sigxx -0.005. Professional 02-08-2026: -0.0075 /
  -0.0099 / -0.009999 for dt 1e-2 / 1e-3 / 1e-4 (the pre-fix GNU gave the
  same digits); the GNU after R6 gives -0.005 for every dt.
- **Root cause**: the Professional appears to keep the legacy
  `tr(new_epe-old_epe)` increment for `-total_linear` (its `-updated_linear`
  gives -0.005 on the same model).
- **Severity**: **LOW** — deliberate divergence; the GNU follows the manual
  (2.2.1: incremental strain = difference of the linear strains) and the
  analytic answer.
- **Proposed fix**: none in the GNU; report to the Professional maintainers.
- **Decision**: divergence accepted by the project owner on 2026-09-29
  (merge `93966ef`); to be investigated further before reporting upstream.
- **Effort**: **cheap**.

### C10. Static pore pressure: the Professional couples it to a post_calcul request

- **Corpus**: Pro `test/other/excavate1` (target node 13 `vely`
  -0.00210987 +- 1e-6).
- **Symptom**: after R7 the GNU gives -0.00253891 (the pre-fix GNU gave
  the target exactly). The Professional gives -0.00253891 too when
  `post_calcul -groundflow_pressure -total_pressure` is added to the input.
- **Root cause**: measured on the Professional 02-08-2026, it applies the
  static pore pressure of `groundflow_phreatic_level_static` to the
  skeleton only when that post_calcul record exists (`-static_pressure`
  does not trigger it). `ground18` without the record: Pro -1745 (no pore
  pressure, wet/dry density still chosen), GNU -1247.5. `excavate1` has
  no such record; its target was produced without the pore pressure.
- **Severity**: **MEDIUM** — one Pro corpus target flips to fail; any
  static-level input without the post_calcul record differs from the
  Professional.
- **Proposed fix**: none by default; the GNU follows the manual (6.574: the
  pore pressures are set to the static pressure; an output request should
  not change the physics). Alternative if exact Pro parity is preferred:
  apply the R7 pore-pressure branch only when that post_calcul record is
  present (one extra check in `groundflow_static_pressure_group()`, for
  the pressure part only).
- **Decision**: divergence accepted by the project owner on 2026-09-29
  (merge `ad3437e`); excavate1 keeps a documented target mismatch.
- **Effort**: **cheap**.

### C11. `control_materi_dynamic` / truss dynamics without `materi_dynamic` diverge from the Professional — RESOLVED 2026-10-11 (F4 T6d)

- **Resolution** (F4 T6d, 2026-10-11): three defects. (1) `so.cc`
  `disp_primary_gate` keyed on the global `materi_dynamic` only and
  excluded `control_materi_dynamic`; now every reader (gate, `solve()`,
  `truss.cc`, `materi.cc`, `elem.cc`) uses `general.cc`
  `materi_dynamic_effective()` (control record at ICONTROL overrules the
  global one, manual 6.141; 1 during the auxiliary f = 1 reassembly), and
  `solve()` no longer PUTs/restores `MATERI_DYNAMIC`. (2) `truss.cc`
  default memory was `-updated`; now `-total_linear` with
  `materi_displacement` (manual 4.12, Pro dynamic3) and `-updated_linear`
  without (Pro probe). (3) the linear branch took `scalar_dabs()` of the
  projected length, wrong for an inverted truss. Results = Pro (node 2
  disx): dynamic3 32.6583 -> 60.44025; only control 0. 6.4477 ->
  0.0958890; global + control 0. 0.0937873 -> 0.0958890; control 0.5
  1.45599 -> 0.932818; no record 1.91225 -> 0.995505; dynamic3
  `-total_linear` -5962.6 -> 60.44025. Tests (test-2026):
  `truss_memory_default_displacement1`, `truss_total_linear_inversion1`,
  `control_materi_dynamic_only1`, `control_materi_dynamic_half1`,
  `truss_dynamic_no_record1`, `truss_memory_default_velocity1`. Residual
  (input errors of the Pro for the memory / materi_displacement
  combinations): C59.

Original report:


- **Corpus**: none; probes = Pro `test/other/dynamic1` (truss, `materi_dynamic
  0.`, control index 10) with extra records (2026-09-29, scratch inputs).
- **Symptom** (node 2 `disx` at t = 100, Pro 02-08-2026 vs GNU):
  `control_materi_dynamic 10 0.` added: Pro 0.0958890 (= dynamic1), GNU
  0.0937873; the same record with index 0: same numbers; only
  `control_materi_dynamic 10 0.` (no global record): Pro 0.0958890, GNU
  6.4477; `control_materi_dynamic 10 0.5`: Pro 0.932818, GNU 0.932764;
  no dynamic record at all (fully implicit): Pro 0.995505, GNU 1.91225
  (with `control_materi_dynamic 10 1.` both give 0.995505).
- **Root cause**: the displacement-primary gate (`so.cc`
  `disp_primary_gate`) keys on the GLOBAL `materi_dynamic` only and sends
  any input with `control_materi_dynamic` at index 0 to the historic path;
  the auxiliary f = 1 reassembly is then overridden by the per-control
  factor read in `truss.cc`/`materi.cc`. The literal index 0 in
  `so.cc:132` is NOT the cause (index 0 and index 10 give identical GNU
  results). The plain implicit truss case (no record) is a separate gap.
- **Severity**: **MEDIUM** — dynamic inputs that use the per-step factor
  or no factor at all differ from the Professional.
- **Proposed fix**: an effective factor per control index (control record
  at ICONTROL, else global) used by the gate, `elem.cc`, `truss.cc` and
  `materi.cc`, with the auxiliary pass protected from the control
  override; then re-verify dynamic1/2/5/8.
- **Effort**: **medium**.
- **Corpus probe** (F2 T2, 2026-10-06): Pro `other/dynamic3` (implicit truss,
  prescribed acceleration, no `materi_dynamic`) is the same gap: node 2
  disx GNU 32.66 vs Pro 60.44 (the target node 1 passes in both).

### C12. bicg "stop" does not abort the GNU run on a singular system

- **Corpus**: none; probe = elasti7 mesh without supports, `solver
  -matrix_iterative_bicg` (2026-09-29, scratch).
- **Symptom**: without `-no` the Professional prints "the BICG solver did
  not converge" and stops without a `.dbs`; the GNU prints the breakdown
  error ("That is not good enough") once per iteration and still finishes
  with rc 0 and a `.dbs` (garbage displacements).
- **Severity**: **LOW** — the error is printed; only the exit status and the
  output differ.
- **Proposed fix**: check the return path of `solve_iterative_bicg()`
  failure in `solve()`/`top.cc` for this breakdown exit.
- **Effort**: **cheap**.

### C13. `control_solver_bicg_stop` at the timestep's index: the Professional ignores it

- **Corpus**: none; probes = elasti7 mesh without supports, `solver
  -matrix_iterative_bicg` (2026-09-29, scratch; same probes as R8).
- **Symptom**: the Professional 02-08-2026 reads only the index-0
  `control_solver_bicg_stop` record, for every timestep; a record with the
  timestep's own non-zero index is ignored (`control_solver_bicg_stop 20
  -no` on a timestep with index 20 still stops). The GNU reads, last one
  wins, `solver_bicg_stop`, the index-0 control record, then the control
  record at the timestep index: 6/9 probes equal the Professional, the 3
  that differ are exactly those with a record at the timestep's non-zero
  index (`20 -no`; `0 -yes, 20 -no`; `0 -no, 20 -yes`).
- **Root cause**: the Professional uses a literal index 0 where the manual
  (6.376) documents a per-timestep index — the same bug class as the
  index-0 `control_materi_elasti_k0` guard fixed in the GNU (R8).
- **Severity**: **LOW** — deliberate divergence; only inputs with a
  non-zero-index `control_solver_bicg_stop` record behave differently, and
  only when the iterative solver fails.
- **Proposed fix**: none in the GNU (the index-0 fallback keeps every
  Professional input that relies on `control_solver_bicg_stop 0 ...`
  unchanged); report to the Professional maintainers.
- **Decision**: divergence accepted by the project owner on 2026-09-29
  (commit pending).
- **Evidence**: `manual-developer/solve_iterative_bicg.md` (3-column
  table), `manual-user/solver.md` (Differences with the Professional
  version).
- **Effort**: **cheap**.

### C14. `element_normal` residual differences with the Professional

- **Symptom**: (1) **RESOLVED** (commit `a8d6452`, 2026-10-01): an
  element whose `element_normal` is zero (interior element, single
  element, middle of a single row) and that uses `_normal_automatic`: the
  Professional aborts with an assertion
  (`!array_null( directions, MDIM*MDIM )`); the GNU now stops with a
  readable input error (it used to apply the full principal-stress
  cut-off). Since 2026-10-03 (commit `b7f3522`) the check is per law: an
  explicit normal of the other law no longer suppresses it
  (`element_normal2d_row_mixed_error1`; Pro: same assertion). (2) quad9/tria6 with curved sides: the
  Professional node normals do not follow the consistent rule (some come
  out axis-aligned); the GNU keeps the rule. (3) the Professional also
  prints `node_normal`, `node_boundary` and `element_boundary`, and
  accepts `global_normal -yes`; the GNU does not. (4) **RESOLVED**
  (commits `d17461e` + `77f9282`, 2026-10-02): the plane direct cut-off
  now follows the Professional rule (separate tension and MC planes,
  tension first with sigy = 0 by default and a principal cap when the
  tension law has no normal, phi_flow return changing only sn and tau),
  reproduced to 4-5 digits on probes c01-c25; the phi_flow factor is a
  deliberate divergence (C15) and the remaining large-strain difference
  of the column is C18 (was C16). Original note: not a normal issue
  but found with it: `mohr_coul_direct 0 0 0` + `_normal_automatic` on the
  upper element of the 2-quad4 column, top sheared by 1 (100
  iterations): both binaries use the normal (0,-1), the shear moves to the
  upper element in both, but the states differ (middle node `disx` GNU
  0.057 vs Pro -0.029; top node `sigxy` 0.036 vs 0, `sigxx` -0.168 vs
  0.043). Pre-fix GNU: `disx` 0.416 (cut-off inert). The plane
  Mohr-Coulomb return (`materi_direct_cutoff`) needs its own comparison;
  no corpus test uses `mohr_coul_direct_normal_automatic`.
- **Evidence**: probes `q_3x3`, `q_row3`, `q9d_2x2`, `t6d_2x2`
  (`manual-developer/element_normal.md`).
- **Proposed fix**: (1) done (input error); (2) none until the
  Professional rule is known; (3) add the records if a user needs them;
  (4) done.
- **Effort**: **cheap**.

### C15. Plane direct cut-off: phi_flow factor M/G instead of the Professional 2 (deliberate divergence)

- **Symptom**: `group_materi_plasti_mohr_coul_direct` with a plane normal
  and `phi_flow != 0`: the Professional 02-08-2026 shifts the normal
  traction by `k (tau - tau_new)` with the constant `k = 2 tan(phi_flow)`;
  the GNU uses `k = (M/G) tan(phi_flow)` (M constrained modulus normal to
  the plane, G shear modulus on it, from the element elastic stiffness).
  Equal at `poisson = 0`. Plane strain, nu 0.25, phi 0.5, c 0, psi 0.3,
  trial sxy 8e-4 / syy -1.2e-3: GNU sxy 7.04154e-4, syy -1.28895e-3;
  Professional 6.92049e-4, -1.26679e-3.
- **Decision**: project owner, 2026-10-02: keep the consistent elastic
  return M/G. Fallback to the Professional constant when no elastic
  stiffness is available (umat).
- **Evidence**: `validation-suite/test-2026/direct_cutoff_plane_phi_flow2.dat`
  (GNU target; the Professional fails it), `_phi_flow1` (nu 0, equal),
  `manual-developer/group_materi_plasti_mohr_coul_direct.md`.
- **Effort**: **cheap** (switch to the constant if convergence is
  preferred).

### C17. `materi_direct_mc` (validation-suite) target stale after C14 — RESOLVED 2026-10-03

- **Resolution**: input corrected to constrain `vely` at all four nodes
  (`bounda_unknown 0 -ra 1 2 3 4 -ra -vely`, the probable intent). The
  target `sigxy` 1 +- 0.3 holds again on the GNU (rc 0), as on the
  Professional with the same constraints. The validation-suite is
  gitignored, so the input change lives in the working tree only.

- **Symptom**: since `d17461e` the plane MC cut-off without
  `tension_direct` caps the principal tension at 0 (Professional rule).
  `materi_direct_mc` constrains `vely` only at nodes 1 and 4
  (`-ra 1 4 -ra` is a list, not a range), so the accumulated compression
  plus the free nodes collapse the element at t = 0.7: `sigxy` 1 for
  t <= 0.6, then 0.129, 0.031, 0.043, 0.057 (target 1 +- 0.3, rc 1). The
  Professional with the same constraints collapses at the same step
  (sigxy 1 then ~0); with all four nodes constrained both binaries keep
  sigxy 1 (sxx = syy -100.87 GNU vs -100.90 Pro with `-updated_linear`),
  and with `-updated` the GNU history follows the Professional
  (1.916..12.06 vs old GNU 2.95..18.69).
- **Proposed fix**: corpus decision — constrain nodes 2 and 3 (the
  probable intent) or recalibrate the target; not done in C14 (the
  validation-suite was outside its edit surface).
- **Effort**: **cheap**.

### C18. Default memory with `materi_displacement`: the Professional behaves as `-total_linear` — RESOLVED 2026-10-07 (F3 T2)

- **Resolution**: `materi_memory_default()` returns `-TOTAL_LINEAR` with
  `materi_displacement` (manual 6.797), adopted together with C8; test
  `materi_memory_default_displacement1` (column disx 0.5 = Pro).
  `mohr_coul_direct7` regression tracked as C41, sfnet `temp5` in B13.

- **Symptom**: after R10 the GNU default is `-updated_linear` for every
  kinematics. With `materi_displacement` and no memory record the
  Professional responds as `-total_linear` (no coordinate update). The
  C14 column (2 quad4, top sheared by 1) still differs: node 3 `disx`
  GNU 0.472 vs Pro 0.5. The GNU updates the coordinates for every memory
  except `-TOTAL_LINEAR` (`elem.cc`, around line 569).
- **Rejected**: mapping the default to `-total_linear` when
  `materi_displacement` is initialized (variant A2) breaks Pro
  `mohr_coul_direct7` and sfnet `temp5`.
- **Proposed fix**: resolve together with the mesh default (C8).
- **Effort**: **large** (tied to C8).

### C19. Explicit `-updated_linear` on a `-follow_material` mesh: reference vs current gradient — RESOLVED 2026-10-07 (F3 T2)

- **Resolution**: the `NODE_START_REFINED` reference-gradient path
  (`c945833`) was removed; tests `materi_memory_updated_linear_current1/2`
  (bar 0.407136 = Pro; column 0.52048 vs Pro 0.52093, in-step geometry
  divergence, see DIVERGENCES-PROFESSIONAL.md).

- **Symptom**: uniaxial stretch to twice the length, explicit
  `group_materi_memory -updated_linear` and `mesh -follow_material`:
  the Professional builds the increment on the current geometry
  (sigxx 0.69565), the GNU on the reference geometry (`NODE_START_REFINED`,
  sigxx 1.0). The default memory already uses the current geometry
  (R10). Measured as a bar in
  `manual-developer/group_materi_memory_updated_linear.md` section 6.
- **Proposed fix**: drop the reference path once the mesh default (C8)
  matches the Professional `-fixed_in_space`.
- **Effort**: **medium**.

### C20. `-updated_jaumann` memory not implemented — RESOLVED 2026-10-08 (F3 T2)

- **Resolution**: `UPDATED_JAUMANN` implemented (linear incremental strain,
  stress rotated with R = I + 0.5(F-F^T), manual 2.2.1); simple shear
  identical to the Professional (tests
  `materi_memory_updated_jaumann_shear1/2`).

- **Symptom**: manual Professional 6.685 lists `-updated_jaumann`; the GNU
  does not accept it.
- **Proposed fix**: add the value with a Jaumann stress rate; verify
  against the Professional in simple shear.
- **Effort**: **medium**.

### C21. `control_materi_updated_apply -yes`: the Professional ignores its manual (deliberate divergence)

- **Symptom**: manual 6.152: `-yes` makes any non-specified memory
  `-updated`. The Professional 02-08-2026 binary ignores it (simple shear
  stays 0 / 0.5 / 0 for displacement, `-follow_material` and
  `-fixed_in_space` kinematics, index 0 and 1). The GNU follows the
  manual: 0.231098 / 0.420739 / -0.228598.
- **Decision**: project owner, 2026-10-04: follow the manual.
- **Evidence**: `validation-suite/test-2026/materi_memory_default_updated_apply1.dat`
  (the Professional fails its targets); no corpus input uses `-yes`.
- **Effort**: **cheap** (drop the `-yes` branch of
  `materi_memory_default()` if the binary behaviour is preferred).

### C22. `hypo12` (Pro `other`) marginal target miss after R10

- **Symptom**: Masin model without a memory record; with the
  `-updated_linear` default `sigyy` misses its target by 1.02e-3 for a
  tolerance of 1e-3 (passed with the old `-updated` default). The file
  belongs to the Professional distribution; its target is not edited.
- **Proposed fix**: none required; revisit with the Masin kernel
  calibration.
- **Effort**: **cheap**.
### C23. `dofid_no` — uninitialized values in the `dof.21` print — RESOLVED 2026-10-05

- **Resolution**: the root cause was NOT an uninitialized `get_new_dbl()`
  buffer but the dof loop of `print_dof()` (`print_hi.cc`): it iterated
  `ipuknwn` over all `nuknwn` value/derivative entries, printed each
  vector/stress block once per entry (shifted into the next dof) and read
  `node_dof[ipuknwn*nder+...]` past the end of the `NODE_DOF` record
  (the next nodes' values and, for the last nodes, heap memory past the
  database array). The loop now runs over the primary unknowns as
  `print_vt.cc` does. `dof.21` 768 lines (20 garbage) -> 48 lines, two
  runs byte-identical; every `dof.*` file is now one line per node per
  component (the documented layout). F0 T6, docs
  `manual-*/print_dof_component_loop.md`. The same loop pattern in the
  FRD/gmsh printers is C24; the uninitialized-buffer audit is C25.


- **Symptom**: `control_print_dof_id -no` (validation-suite `dofid_no`)
  writes 20 of 768 lines of `dof.21` with heap garbage (`6.9e-310`,
  different on every run) at node (1,1). Found by the F0 T4 old/new A/B
  (2026-10-05); present in both binaries, unrelated to the
  `generate.cc` `element_dof` fix. The run and its file check (3
  columns) pass.
- **Proposed fix**: find the uninitialized buffer in the `print_dof()`
  path for the `-no` id layout (same pattern as the generated
  `element_dof`: a `get_new_dbl()` without `array_set()`).
- **Severity**: **LOW**. **Effort**: **cheap**.

### C24. `print_frd` / `print_gmsh` — repeated blocks and out-of-record reads — RESOLVED 2026-10-05

- **Resolution** (F1 T1): `print_frd()` and `print_gmsh()` now loop over
  the primary unknowns like `print_vt.cc` / `print_dof()` (type, label
  and dof type at `iuknwn = ipuknwn*nder`, step `nval`). FRD 2D vector
  blocks write 0 as third component (it read the next dof); gmsh field
  names are built from the dof type (`node_materi_velocity_0`, Pro base
  names). Test `test-2026/print_frd_gmsh_dof_loop1` (quad4, `derivatives`):
  68 -> 5 FRD blocks, 324 -> 17 gmsh `node_` fields, values equal to the
  first block of each dof of the Professional 02-08-2026; checked in
  `build_safe.sh`. The existing `freq_timeint`/`freq_timestep` FRD checks
  only read the `100CL` time and are unaffected. Remaining layout
  differences: C27. Docs `manual-*/print_dof_component_loop.md`.


- **Symptom**: `print_fr.cc` (`print_frd`, result blocks, line ~229) and
  `print_gm.cc` (`print_gmsh`, node data, line ~442) use the same loop as
  the old `print_dof()` (C23): `for ( ipuknwn=0; ipuknwn<nuknwn_;
  ipuknwn++ )` with `dof_scal_vec_mat[ipuknwn]` and `base_indx =
  ipuknwn*nder_`. Each vector/stress result block is written once per
  stored value/derivative entry, and the indices run past the end of the
  `NODE_DOF` record (garbage for the last nodes, worse with
  `derivatives`).
- **Proposed fix**: the primary-unknown loop of `print_vt.cc` /
  `print_dof()` (`npuknwn = nuknwn/nder`, step `nval`, type from
  `dof_scal_vec_mat[iuknwn]`). Changes the number of FRD/gmsh blocks: the
  `freq_*` and gmsh checks of `build_safe.sh` need review.
- **Severity**: **LOW** (output files only). **Effort**: **cheap**.

### C25. Uninitialized-buffer audit follow-ups (F0 T6)

- **Audit**: every `get_new_dbl(` / `new double[` site in `*.cc` was
  checked for the `a32ea74` pattern (buffer partly assigned, then PUT or
  printed): none confirmed. Uncertain or related items, not fixed:
  1. `calcul_force.cc` ~2824: `pef_vel_prev` is reallocated larger
     without zeroing (only the first allocation is zeroed); a node
     inactive at refresh and active later would read garbage `vprev`
     into the `post_element_force` inertia share. Fix:
     `array_set( pef_vel_prev, 0., need )` after the reallocation.
  2. `print_te.cc` ~151 and `print_pl.cc` ~436: `NODE_EIGEN` read with
     `GET_IF_EXISTS` into a reused buffer; a node without the record
     prints the previous node's values (so.cc writes it for every active
     node, so probably unreachable).
  3. `split.cc` ~55-68: the zeroed `ddum` is also the output array of
     `geometry()`, so the "zero" fallback passed to `create_node` may hold
     geometry results (contamination, not uninitialized memory).
  4. Side findings: `print_beam_force_moment.cc` ~282 passes a
     one-element `ddum` as the 3-element work array of `array_distance`
     (stack overflow); `interface.cc` ~1044 `free()`s `new[]` memory;
     `mpc.cc` ~1050 reads `MPC_NODE_FACTOR` into an `nmaster` buffer
     before clamping its length; `nonloc.cc` ~326 PUTs a stale
     `ielem_dof` when `ielem_nonloc==0`.
- **Severity**: **LOW** (edge cases, none observed in the corpora).
  **Effort**: **cheap** each.

### C26. `control_print_dof -separate_index` file layout differs from the Professional

- **Symptom**: `dofid_no` (converted to the Pro dialect) with
  `control_print_dof 21 -separate_index`: the Professional 02-08-2026
  writes one file per dof (`disx.21`, `disy.21`, `sigxx.21`, ...,
  `coord_0.21`, `coord_1.21`), each line `x y value` in `%e` with 12
  digits; the GNU writes one `dof.21` with the components one after
  another (`x y value`, short `%g`). The values per component agree
  (e.g. `disy` 0.002 at the top nodes) after `f22ba98` (C23); before it
  the GNU file also repeated blocks and read past the record.
- **Proposed fix**: follow the Professional layout (one file per dof name,
  Pro number format) or document the GNU layout as a divergence;
  check `control_print_dof_id` interplay and the build_safe `dof.*`
  checks.
- **Severity**: **LOW**. **Effort**: **cheap/medium**.

### C27. FRD / gmsh file layout differs from the Professional

- **Symptom** (found by F1 T1, `print_frd_gmsh_dof_loop1` converted and
  run on the Professional 02-08-2026): the dof values agree, the layout
  does not. FRD: the Professional writes a `1PSTEP` + `100CL` header per
  result block, vector blocks with 4 components (`V1 V2 V3 ALL`, all
  with `icind1 = 1`), the names `VELOCITY` (GNU `VELO`), `S11..S13`
  (GNU `SXX..SZX`), `T11..T13` (GNU `EXX..EZX`), extra non-dof blocks
  (`BOUNDA_X/Y/D`, `BOUNDARY`, `REACTION`), its own header lines
  (`1UDATE`, `1UTIME`, `1UPGM CalculiX`), and with `derivatives` one
  block per value AND derivative entry, whose values are shifted into
  other dofs (the C24 defect class, in the Professional). gmsh: the
  Professional writes vectors/tensors as ONE 3/9-component field
  (`node_materi_velocity`, `node_materi_stress`), element fields from
  the integration points (`element_materi_stress`), plus
  `node_bounda_dof_*`, `node_boundary`, `*_rhside`, `element_group`,
  `element_flux_temperature`; the GNU writes one scalar field per
  component (`node_materi_velocity_0`) and nodal averages per element.
- **Proposed fix**: decide per format which Professional conventions to
  follow (post-processors: CalculiX cgx, PrePoMax, FreeCAD, gmsh); do NOT
  copy the derivative-block defect.
- **Severity**: **LOW** (output files only). **Effort**: **medium**.

### C28. `node_force` applied with the opposite sign (and doubled in a 1D probe) — RESOLVED 2026-10-06 (F1 T6)

- **Symptom** (found by F1 T3, K04 `force18` probes on the Professional
  02-08-2026): Pro `force18` with `node_force 2 1.0` instead of the
  `force_point`: Pro node 2 `velx` +1, GNU −2. `dof.cc`
  `parallel_new_dof_before()` does `node_rhside -= node_force` (the
  `force_point` had the same sign error, fixed in F1 T3), and the
  magnitude is doubled in this 1D bar2 model (cause not traced).
  `delete.cc` (mesh_cut node forces) negates the stored value to
  compensate the sign (documented there).
- **Proposed fix**: add the force (`+=`), drop the `delete.cc`
  negation, trace the factor 2, A/B the corpus inputs with `node_force`.
- **Resolution** (F1 T6): `dof.cc` adds the force; `delete.cc` stores the
  mesh_cut substitute without the negation (mesh_cut_1 `node_force 6` +1
  = Pro). The factor 2 is not a double application: the reversed load
  compressed the unit bar by exactly its length (singular zero-length
  element; `node_force 2 -1.0` still gives GNU -2 / Pro -1, `-0.5` both
  -0.5). Probes 1D/2D/3D = Pro. Corpus A/B: only `node_force_inertia`
  changes (stale target, B12, retargeted to -20.11).
- **Severity**: **MEDIUM** (wrong load direction). **Effort**: **cheap/medium**.

### C29. `group_condif_conductivity` with one value per direction not accepted

- **Symptom** (F1 T3, Pro test `condif11`): `group_condif_conductivity 0
  0.1 0.1` is a PARSE error in the GNU ("I do not know : 0.1"); the
  Professional accepts it (anisotropic conductivity, as
  `group_groundflow_permeability`).
- **Proposed fix**: accept ndim values, isotropic when one value.
- **Severity**: **LOW**. **Effort**: **cheap/medium**.

### C30. `mohr_coul_direct4` (Pro `other`): the Professional does not enforce the interface Mohr-Coulomb strength (deliberate divergence)

- **Symptom** (F1 T3, audit K04): interface shear test (phi 0.2 rad,
  c 1, normal stress -1, top displaced 1 at t = 1 and held to t = 100).
  Manual Professional 6.632: maximum friction c + |Fn| tan(phi) = 1.2027.
  The GNU returns to it at once (sigxy 1.2027 from t = 1, `disx(3)`
  0.3007, status plastic). The Professional 02-08-2026 stays at the
  elastic 2.0 with status `-plastic` (t = 1: 1.99994, t = 5: 1.9911,
  t = 20: 1.8283, t = 50: 1.0953) and ends below the strength at
  t = 100 (0.40976, status `-elastic`, target `disx` 0.10244).
- **Root cause**: Professional behaviour (no immediate return, slow
  drift), PRO-suspect; the GNU follows the manual.
- **Decision** (project owner, 2026-10-06): deliberate divergence. The GNU
  keeps the manual 6.632 immediate return to the strength; the
  Professional's slow drift looks like an artefact of its own integration
  and is not imitated. Recorded in "Differences with the Professional
  version" of `manual-user/` and `manual-developer/group_materi_plasti_mohr_coul_direct.md`.
- **Severity**: **LOW** — deliberate divergence; the GNU follows the manual.
  **Effort**: —.

### C31. `tet10` integration rule is only degree 1 (Pro `tet10_2` under gravity) — RESOLVED 2026-10-09 (F4 T4)

- **Resolution (F4 T4, option (c), project owner 2026-10-09)**: the
  Professional's Keast 5-point rule (`polynom.cc`), its linear
  least-squares nodal recovery `1,L1,L2,L3` (`miscel.cc`) and a positive
  HRZ lumping of the principal-dof inertia/capacity (`general.cc`,
  corner V/36, mid-edge 4V/27; deliberate divergence TET10-HRZ: the
  Professional row sum gives the corner -V/20 and its dynamics fail).
  `tet10_2` sigzz -1.00000000515 (Pro -1); tests
  `test-2026/tet10_body_load1`, `tet10_mass_hrz1`,
  `tet10_dynamic_gravity1`. Found on the way: the Professional
  `post_point_dof` inside a tet10 interpolates linearly from the corners
  (DIVERGENCES TET10-PPOINT, pending decision).


- **Symptom** (F1 T3, audit K04): `tet10_2` (hex8 -> h -> p -> split into
  48 tet10, gravity): GNU post-point sigzz -0.894 (nodal sigzz at z = 0
  between -0.94 and -0.88, top disz -0.48..-0.56), Professional exact -1
  / -0.5. A single hex27 (no split) is exact in the GNU; `tet10_1`
  (prescribed strain, no body force) passes.
- **Root cause**: `polynom.cc` TET10 branch: nodal quadrature with weight
  ~0 (EPS_TRIANGLE) at the corners and 1/6 at the 6 mid-edge nodes. On a
  tetrahedron this rule is exact only for linear integrands (the
  triangle analog, mid-edge 1/3, is degree 2): the quadratic stiffness
  and the consistent body load (corner share -V/20, mid-edge V/5) are
  wrong.
- **Experiment** (scratch build, not committed): the degree-2 nodal rule
  (corners -1/20, mid-edge 1/5) makes `tet10_2` exact (-1.0000000058),
  `tet10_1` stays exact.
- **Options**: (a) that nodal rule: exact, but NEGATIVE corner weights,
  so the lumped inertia/capacity (general.cc, h = delta at nodal points)
  becomes negative at the corners; (b) a 4/5-point Gauss rule: breaks
  the npoint == nnol assumption of the element_dof storage; (c) (a) for
  the stiffness/loads plus a separate positive lumping (HRZ) for the
  inertia. Needs a design decision and an A/B on the sfnet tet10 inputs
  (`tet10`, `refine4`, `ho_othr1/2`, `examp*`).
- **Decision** (project owner, 2026-10-06): moved to **F4** (design +
  corpus A/B); not part of F1.
- **Severity**: **MEDIUM** (wrong quadratic-tet results under body
  loads/bending). **Effort**: **medium**.

### C32. K04 targets that are not GNU defects (`axisym9`, `distri2`) — CLOSED 2026-10-09 (F4 T4)

- **Closure (F4 T4)**: neither target is a GNU defect and nothing is left
  to do here.
  - `axisym9`: the row was stale. The GNU already gave post_point sigxx
    0.00159154943 = Professional (0.001591549431; target 0.0015915 +/-
    1e-6) before F4 T4; the 0.00158902 of F1 T3 predates the Gauss
    alignment of the recovery. The forced 1-point axisymmetric rule
    (audit K20) it was attributed to was removed in F4 T4 and the value
    stays 0.00159154943.
  - `distri2` (post_point temp 0.998933 vs 0.999774 +/- 1e-4):
    `control_distribute -normal` random field; both binaries draw a valid
    N(1, 2e-3) (GNU mean 0.99995 / std 0.00203, Pro 0.99997 / 0.00201)
    from different random generators. The target encodes the
    Professional random sequence; no fix (DIVERGENCES DISTRI, class C).

### C33. `_visco` direct cut-offs: one-step-lagged relaxation time step and no visco on the principal tension cap — RESOLVED 2026-10-09 (F4 T3)

- **Resolution**: `set_stress()` reads `DTIME` at `VERSION_NEW`; the
  principal tension cap applies the tension `_visco` factor; in the full
  principal mode each law relaxes with its own `_visco` record. 12 direct
  probes and 2 incremental viscoplastic probes equal to the Professional
  02-08-2026; tests `materi_direct_visco_first_step1/2`,
  `materi_direct_visco_separate1`, `materi_plasti_visco_power_dt1`.

- **Symptom** (F1 T4, K13 probes on the Professional 02-08-2026, single
  quad4, E 1, nu 0.25, tension_direct 1e-3, `_visco` tm 2, eps_xx 2e-3 in
  one step dt 1): Professional sigxx 1.84914e-3 (= 2.4e-3 - (1 -
  exp(-0.5)) 1.4e-3), GNU 1e-3 (full cut). `compression_direct_visco`
  (same data, compression): Professional -1.84914e-3, GNU -2.4e-3 (no
  cut at all); `mohr_coul_direct_visco` full and plane mode: GNU no cut
  in the first step (sigxy 2.4e-3, Professional 1.7575e-3 / 1.7996e-3).
- **Root cause**: `stress.cc` `set_stress()` reads `DTIME` from
  `VERSION_NORMAL` (the previous step; 0 in the first step, so the
  factor 1 - exp(-dt/tm) is 0), every other module reads `VERSION_NEW`;
  and the principal tension cap (`direct_cap_principal_values`, full mode
  and plane mode without a tension normal) ignores the `_visco` factor.
  `visc_pl4` (1D, 200 steps) still ends at the target 1 in both.
- **Proposed fix**: `VERSION_NEW` in `set_stress()` (also used by the
  incremental viscoplastic laws: `lambda*dtime`, needs its own corpus
  A/B) and the visco factor on the principal cap.
- **Decision** (project owner, 2026-10-06): scheduled for **F4** (the
  `VERSION_NEW` change also touches the incremental viscoplastic laws and
  needs a corpus A/B); not part of F1.
- **Severity**: **MEDIUM** (`_visco` direct laws relax with the wrong
  time step). **Effort**: **cheap/medium**.

### C34. Direct cut-off shear plastic strain: the Professional writes twice the tensor value (deliberate divergence)

- **Symptom** (F1 T4, K13): since K13 the GNU stores the direct cut as
  plastic strain C^-1 (sigma_trial - sigma). Normal components and the
  yield value equal the Professional; the shear components differ by a
  factor 2: full MC cut `direct_plastic_strain_mc1` GNU eppxy
  2.04114892e-3 = dsig_xy / 2G (tensor), Professional 4.08229785e-3 =
  dsig_xy / G. The Professional's own total strain (eptxy) and its
  incremental-law plastic strain (von Mises shear probe: eppxy
  2.27831e-3 in both) are tensor components, and its epexy = eptxy -
  eppxy then violates sigma_xy = 2G epe_xy.
- **Decision** (project owner, 2026-10-06): **deliberate divergence**.
  The Professional writes the shear plastic strain of a direct cut as
  dsig / G (twice the tensor component) while its total strain uses
  tensor components; the GNU keeps the tensor value so that
  sigma = C (eps - eps_p) holds. Documented in "Differences with the
  Professional version" of `direct_plastic_strain.md` (user and
  developer manuals) and in the SEGUIMIENTO "Descarte de features"
  table.
- **Severity**: **LOW** (output only). **Effort**: —.

### C35. One-equilibrium-iteration lag of the whole material state (audit K02) — RESOLVED 2026-10-10 (F4 T5), gated families open (C52, C53)

- **Symptom** (F1 T5, Professional 02-08-2026, sfnet `spring4` converted:
  rigid translation velx 1, `materi_strain_total`, one equilibrium
  iteration): GNU node and integration-point eptxx = sigxx = -0.01,
  element_dof velx 0.5; Professional 0 / 0 / velx 1. With 3 iterations
  both give 0.
- **Root cause**: the constitutive evaluation (and the stored
  integration-point state ELEMENT_DOF) of an iteration uses the
  PRE-solve velocity; the hoisted post-solve pass of `top.cc` only
  re-assembles the nodal STRESS and only with `stress_recovery_absolute`
  (off with `materi_strain_total`), so the strain/stress state of the
  last iteration lags the solved velocity. The Professional evaluates the
  material on the post-solve velocity (ip velx 1).
- **Why parked**: a fix is a post-solve material pass that also updates
  the integration-point history (ELEMENT_DOF) and the incremental
  non-principal dofs (general() denominator/inertia), i.e. a change of
  the iteration scheme for every model with few iterations; needs a
  design decision and a full corpus A/B (the audit counted 6 tests that
  vanish at 5 iterations and 10 that improve).
- **Decision** (project owner, 2026-10-06): scheduled for **F4** (design +
  full corpus A/B); not part of F1.
- **Severity**: **MEDIUM** (one-iteration models report a lagged state).
  **Effort**: **medium/large**.
- **Resolution** (F4 T5, project owner 2026-10-10): post-solve material
  state pass as the default (`top.cc`, after the equilibrium loop:
  `iteration_start(); element_loop(); parallel_new_dof_diagonal;
  mpc_node_apply()`; `material_state_post_pass` makes `elem.cc` skip the
  matrix store). It replaces the stress-only sync except for the gated
  families: inertia, slide, contact (temporary divergence, DIVERGENCES
  class B `POSTPASS-GATE`). Verified: probes = Professional (`spring4`
  eptxx 0, plastic bars eppxx 0.01, stretch velx 0.5 every step); matrix
  vs no-matrix byte-identical; corpus A/B: see SEGUIMIENTO F4 T5 row.
  Exposed and fixed first: C51 (hardsoil kappa read). Docs:
  `manual-*/post_solve_material_state_pass.md`.
- **Widening criterion**: a gated family is enabled only after the pass
  is extended to its physics (`slide()` / `parallel_contact` called
  consistently, consistent inertia update) and the probes + full corpus
  A/B on omen show no test moving away from the Professional.

### C36. `group_materi_elasti_young_power` next to `group_materi_elasti_young`: the Professional adds both stiffnesses — DECIDED 2026-10-06: deliberate divergence, option (a) (audit K06)

- **Symptom** (F2 T1, Professional 02-08-2026, validation-suite `myoung6*`
  converted): `myoung6_e2` (law clamped to E2 = 3000, `young` 1000) Pro
  sigyy -5.3846 (stiffness 4000), GNU -4.0385 (3000); `myoung6_e3` (law
  clamped to 800) Pro 1800, GNU 800; `myoung6` Pro -3.3485, GNU -1.6657.
- **Measured**: without the `group_materi_elasti_young` record the Pro
  equals the GNU on all three (-4.0385, -1.0769, -1.6742 vs -1.6657).
  With `materi_history_variable 1` the Pro reports the "calculated young"
  (manual 6.661) as 3000 on `myoung6_e2` while its stress needs 4000: the
  Professional builds C(young) + C(E_law), the same accumulation the GNU
  2014 code had (removed in Sprint 10).
- **Manual**: theory 2.2.2 "the young's stiffness modulus is made a
  function of the average stress state: E = E0 + E1 (p/p1)^alpha" and the
  history variable of 6.661 both say the law IS the modulus; the manual
  is silent on combining it with a separate `group_materi_elasti_young`.
- **Options**: (a) keep the GNU (law = modulus, manual reading; the
  `myoung6*` targets stay; documented divergence from the binary);
  (b) follow the binary (add C(young) when both records exist; retarget
  `myoung6`, `myoung6_e2`, `myoung6_e3`; `myoung6_apply` with
  `materi_elasti_young_power_apply -no` must then be re-measured on Pro).
  Only the validation-suite `myoung6*` family and Pro `validation_14_mesh`
  use the record. `group_materi_elasti_poisson_power` is NOT affected:
  there the Pro replaces the poisson ratio (fixed in F2 T1).
- **Decision (project owner, 2026-10-06)**: option (a), deliberate
  divergence. The GNU follows the manual (theory 2.2.2, 6.661: the power
  law IS the modulus); the code and the `myoung6*` targets stay. Recorded
  in SEGUIMIENTO "Descarte de features" and in "Differences with the
  Professional version" of `manual-user/` and
  `manual-developer/group_materi_elasti_young_power.md`. No code change.
- **Severity**: **LOW** (input with both records only). **Effort**:
  **cheap** either way.

### C38. `control_support_edge_normal_stiffness_freeze` is parse-only

- **Manual** 6.381: with `-yes` the stiffness FORCES of
  `support_edge_normal` are frozen at their present value (gravity with
  support stiffness first, then vibration with damping only).
- **GNU**: parsed and accepted, no effect (the earlier note "meaningful
  only with plasticity" misread the manual).
- **Measured** (F2 T2, Pro `support_edge_normal_3`): Professional disy
  5.00999e-4 with the record, 5.00998e-4 without; GNU 5.00998e-4 (after the
  damping-tangent fix). The target (tolerance 1e-6) passes either way.
- **Proposed fix**: at the control index with the record, store the
  per-node stiffness force `k u` (normal and tangential) and use it
  instead of `k u` (and drop `k*dt` from the matrix) while frozen.
- **Severity**: **LOW**. **Effort**: **cheap/medium**.

### C37. 3D beam: the GNU beam is planar, the Professional beam is a full 3D frame element (audit K17 `trubea2/3`) — PARKED, design

- **Symptom** (F2 T2, Pro `other/trubea2`, `trubea3`): a truss-beam along x
  loaded in z (x-z plane) or in the y-z plane: GNU bicg breakdown, disz
  2e5; Professional 02-08-2026 disz 0.00333 and roty/rotx reactions
  exact. `trubea1` (x-y plane) and `trubea4` agree.
- **Root cause**: `beam_3d()` (`beam.cc`) maps a 2D beam into ONE plane
  given by the GNU-only `group_beam_plane` (default -x -y); out-of-plane
  bending has no stiffness. The Professional (manual 6.590-6.592) has a
  3D beam with local axes: `group_beam_inertia Iyy Izz J` (two bending
  planes + torsion), `group_beam_direction_z` (default: an arbitrary
  normal in 3D), output `element_beam_direction`.
- **Why parked**: a new 12-dof 3D frame element (local triad, torsion,
  both bending planes, -updated large rotations, beam force/moment output
  in local axes, `print_beam_force_moment`), not an isolated fix.
- **Severity**: **MEDIUM** (any 3D beam out of the x-y plane). **Effort**:
  **large**.

### C39. K17 residuals: inactive-element nodes, `update2`, `truss13` (audit K17)

- **`truss8`** (Pro `other`): a truss born at t = 1 (`group_time`), node 1
  driven by a constant velocity 0.1 from t = 0. Node 2 disx 0.1 in both
  (target); node 1 disx GNU 0.2, Professional 0.1: the Professional does
  not move a node that has no active element (manual 6.769 is silent).
  **RESOLVED 2026-10-06** (project decision: follow the Professional;
  F2 T2c, hash in `odd/tasks/f2-medium-bugs.md`): `bounda.cc`
  `bounda_node_dormant()` skips `bounda_dof`/`bounda_unknown` on a node
  whose elements are all inactive by `group_time`; truss8 node 1 disx
  0.1 = Pro at every step; test `group_time_bounda_dormant1`. **Still
  open (LOW)**: after the DEATH of a group the Professional brings its
  nodes to rest (velocity 0, also the free ones; probe `group_time 0 0.
  1.`), the GNU keeps their old velocity (no equation); and only
  `group_time` makes a node dormant (mesh_activate_gravity method 1 and
  mesh_delete not measured).
- **`dynamic3`** (Pro `other`): implicit truss dynamics without
  `materi_dynamic`; the target node passes, node 2 disx GNU 32.66 vs Pro
  60.44. Same gap as C11 ("plain implicit truss case"); tracked there.
- **`update2`** (Pro `other`): NOT a truss case (misclassified in K17): a
  -updated quad cantilever 100 x 10, E = 1e7, edge force 1e6 with
  `force_edge_normal_factor`, p-refinement. GNU "solution seems to
  diverge" (veliy -86434), Professional veliy 12.243 (target 12.24).
  Without p-refinement or with -updated_linear both stop. Large-rotation
  robustness of the -updated material, related to K19; not investigated
  further.
- **`truss13`** (Pro `other`): `control_mesh_truss_distribute_mpc` is not
  implemented (SEGUIMIENTO "PENDIENTE"); the truss across an interface
  carries no force in the GNU (Pro 5.0).
- **Zero-force straight truss**: the first iteration of a straight truss
  under a transverse load is still singular (no force, so no geometric
  stiffness); the GNU prints "iterative solver broke down" and continues
  (C12), then converges (Pro `truss2/3/14` equal after F2 T2).
- **Severity**: **LOW**. **Effort**: **cheap** (truss8 decision) to
  **large** (update2).

### C40. Element-target contact residuals (audit K05, after F2 T3b)

- **Fixed in F2 T3b (2026-10-07)**: with `contact_target_element_group`
  (manual 6.105) the target sides need `node_boundary` (6.875); the GNU
  never generated it (only `control_mesh_macro`, `adjust`, `delete` and
  hand-written input did), so the Professional inputs `contac1`, `contac4`
  and `axisym3` had no contact. `area.cc` `contact_node_boundary_set()`
  now generates it (nodes of unshared sides) when the group record is
  given and no node_boundary exists; contac1/contac4 = Pro exactly,
  axisym3 velx -9.99981e-4 (= Pro), 2D quad4 probe = Pro.
- **Open**: (1) the Professional refuses a target element whose nodes are
  all on the mesh boundary ("Error, a target element has all nodes on the
  boundary of the mesh", manual 6.105 "more than 1 layer"); the GNU runs
  it and gives a non-physical result (single-column probe: the far nodes
  move). (2) `parallel_contact` accepts a target side when
  `number_of_side_boundary_nodes==ndim`: a hex8 side has 4 nodes, so hex8
  targets are never detected (code reading, not measured on the Pro). (3) The
  Professional writes `node_boundary`/`element_boundary` in every run
  (also without `mesh_boundary`, against its manual 6.815 default -no);
  the GNU generates node_boundary only for element-target contact and has
  no `element_boundary` record.
- **Severity**: **LOW**. **Effort**: **cheap** (1, 2) to **medium** (3).
---

### C41. `mohr_coul_direct7` (Pro `other`): `-total_linear` with a velocity `mpc_node_number` tie (F3, 2026-10-07) — tie RESOLVED 2026-10-09 (F4 T2); residual open

- **Resolution (F4 T2)**: root cause `mpc.cc` `mpc_node_apply()`: the
  post-solve re-sync set the slave VELOCITY but not its integrated
  displacement (`parallel_new_dof_diagonal` integrates `dis = dis_old +
  vel*dt` before the re-sync), so the slave displacement lagged one solve;
  with `-total_linear` the tied nodes 3/4 ended with disy -0.001808 /
  -0.001692. The slave displacement is now integrated with the slave
  velocity. Target met: sigxy 0.0998334 at every node (= Professional);
  tied nodes equal. Test `test-2026/mohr_coul_direct_mpc_total1`.
- **Residual (open, LOW)**: disy -0.00175 vs Professional -0.0018458 and
  sigzz -1.0084 vs -0.90785 (sigxx, sigyy, sigxy equal). Not the tie (both
  programs keep the tied nodes equal); belongs to the `mohr_coul_direct`
  cut under `-total_linear`, possibly through the stored direct plastic
  strain (C34), not verified. DIVERGENCES C41.

- **Corpus**: Pro `other` `mohr_coul_direct7` (exit-code pass before F3,
  fail after). The input has `materi_displacement` and no memory record,
  so since F3 (C18) the default memory is `-total_linear`, as in the
  Professional. With `-total_linear` and the velocity `mpc_node_number`
  tie of the input the GNU misses the target: node 1 sigxy GNU 0.102732
  vs target 0.0998334 (Professional meets it; GNU before F3, default
  `-updated_linear` with coordinate update: 0.0998328).
- **Root cause**: not established; a separate GNU defect in the
  `-total_linear` path with a velocity tie (the same input with an
  explicit `-total_linear` already failed before F3).
- **Severity**: **LOW**. **Effort**: **medium** (investigation).
### C42. quad4 nodal stress recovery: natural-coordinate bilinear interpolation (GNU) vs physical-coordinate linear least squares (Professional) (F4 T4, 2026-10-09)

- **Symptom** (F4 T4 investigation, after the K20 fix): on rectangles in
  plane analyses both programs agree; on distorted quad4 and on
  axisymmetric rings the nodal stresses differ although the
  displacements and the integration-point stresses are equal.
  Distorted 4x4 axisymmetric patch `hg4d`: nodal stresses 3-7 % apart;
  Akin ring (`axisym2`, test `test-2026/axisym_quad4_ring1`): inner sigzz
  GNU 10760.09, Professional 10759.22 (disx equal to 1e-12).
- **Root cause**: `sri_stress_recovery_weight()` (`miscel.cc`) recovers
  the quad4 by the bilinear `1,xi,eta,xi*eta` interpolation of the 2x2
  Gauss values in the natural coordinates (identified on a rectangle,
  2026-09-13). The Professional fits `1,x,y` by uniform-weight least
  squares in the PHYSICAL coordinates of the points and evaluates it at
  the physical node (the two coincide on a rectangle up to the xy term).
- **Experiment** (variant v4 of the F4 T4 investigation, scratch build):
  the physical-coordinate fit reproduced the Professional `axisym2`
  sigzz (10759.2155) and moved the distorted `hg4d` nodal stresses
  towards the Professional, but changed the exit code of `cd_method`
  (validation-suite, pass -> fail), so it was not adopted with K20.
- **Severity**: **LOW** (post-processing of the nodal stress field only;
  the solution is unchanged). **Effort**: **medium** (the per-element
  physical fit, plus the `cd_method` analysis).

### C43. Memory hogs of the default Bi-CG solver (F4 T7, 2026-10-09) — memory RESOLVED; solver speed RESOLVED (F4 T9, 2026-10-10)

- **Symptom**: five Professional `other` inputs (`force11`, `large1`,
  `ground19_water_under_dam`, `mesh_activate1`,
  `foundation_consolidation_settlement`) need more than 4 GB with the GNU
  (bad_alloc under the cap, up to 11.9 GB RSS under a 12 GB cap) while
  the Professional runs them in 0.07-1.3 GB.
- **Root cause A** (`so.cc` `solve()`): since `5b7bb8b` every solve
  allocated and zero-filled a dense `n x n` matrix `mat_dense` that was
  never read (its `dgesv_` reader was a dead branch). **Resolved** in
  F4 T7 A: removed; no result change.
- **Root cause B**: the band LU workspace of the Bi-CG retry was allocated
  and filled on every solve with the raw node-number band. **Resolved** in
  F4 T7 B: built lazily in `solve_band_retry()` after a failure, exact
  band, Cuthill-McKee ordering when narrower; round-off changes in the
  retried solves only.
- **Root cause C** (`mesh_activate1`): `control_mesh_macro` on the
  `control_timestep` index added a mesh copy per step. **Resolved** in
  F4 T7 C (C44).
- **Measured** (omen, 12 GB cap, old -> new, Professional): `force11`
  bad_alloc 7.8 GB -> 0.27 GB (Pro 0.27 GB); `large1` bad_alloc -> rc 0,
  0.65 GB, 58 s (Pro 1.3 GB, 9 s); `ground19` 5.2 -> 0.12 GB (Pro 0.17);
  `mesh_activate1` 11.9 GB bad_alloc -> 0.035 GB, 29 s (Pro 0.07 GB, 26 s);
  `foundation_consolidation_settlement` bad_alloc -> 1.5 GB, timeout 600 s
  (Pro timeout too). Table: `manual-developer/solver_memory.md`.
- **Residual (solver speed)** after F4 T7: `force11` (one 3D solve)
  did not converge with the Bi-CG within 600 s (Pro PARDISO 13 s);
  `foundation_consolidation_settlement` breaks the Bi-CG down on every
  solve and runs the band LU retry (1.5 GB).
- **F4 T9 (2026-10-10, owner decisions of the same day)**. Findings: the
  Professional default solver is PARDISO, and the Professional Bi-CG does
  not converge on `force11` either; the GNU `-matrix_pardiso` fell back to
  the Bi-CG; the linked SuperLU 6.0.1 route crashed (use-after-free in
  `so_suplu.c`); the Bi-CG products were element by element; the nodal
  stress recovery re-solved a 4x4 least-squares system 43.7 M times on
  `large1`. Changes: `681dc83` SuperLU use-after-free and new[]/free
  mismatch; `789f8bf` the same db-pointer-into-PUT use-after-free in
  `interface.cc` `interface_convert()`; `acfac40` recovery factor table
  (bit-identical); `3349515` `-matrix_pardiso` -> SuperLU with Bi-CG
  fallback and memory guard (GNU default stays the Bi-CG: SuperLU does
  not scale in 3D, `large1` 30^3 > 4 GB); `b7009dd` CSR operator for the
  CG/Bi-CG products (~5x per iteration in 3D). C12/C13 `bicg_stop`
  semantics unchanged. `force11` with its own `solver -matrix_pardiso` now
  finishes in 15 s (0.41 GB; target miss C47); `large1` 60 -> 24 s;
  corpus -32 % (tests finished on both).
  Measurements: `manual-developer/solver_memory.md` (F4 T9 table).
  Remaining: `foundation_consolidation_settlement` (default Bi-CG,
  breakdown + band retry on every solve; the Professional times out too).

### C44. One-shot mesh records re-applied on every time step (F4 T7, 2026-10-09) — RESOLVED

- **Symptom**: `control_mesh_macro`, `_copy`, `_mirror`, `_switch` on the
  same index as a multi-step `control_timestep` were applied on every
  step (`mesh_activate1` 202 -> 4040 nodes in 20 steps; copy doubled the
  mesh per step; switch swapped the axes back).
- **Fix**: `top.cc` `step_start_control_entry` / `first_of_control`: run
  them on the first step of the index (Professional probes: once per
  index). Tests `mesh_macro_timestep1`, `mesh_switch_timestep1`.

### C45. `control_mesh_copy` / `control_mesh_mirror` after the first time step stop with a `node_dof_start_refined` error (F4 T8, 2026-10-10) — RESOLVED

- **Symptom** (F4 T7 probes): a copy or mirror on an index after a
  `control_timestep` index stopped with "Error detected for data item :
  node_dof_start_refined / record : 8" (rc 1); the Professional runs it.
  Also on index 1 without an earlier time step.
- **Root cause**: `mesh.cc` `mesh_copy()` / `mesh_mirror()` copied only
  `NODE_DOF` to the new nodes; `renumber()` reads `NODE_DOF_START_REFINED`
  of every node. With that fixed, `create_element()` reset the
  integration-point state of the duplicates (1 % off the Professional on
  a loaded block).
- **Fix** `0ca144c`: `mesh_duplicate_node_data()` (all node-class items,
  buffered PUT) and `mesh_duplicate_element_state()` (`ELEMENT_DOF`,
  `ELEMENT_DOF_INITIALISED`); element numbers `old + max + 1` (Pro).
  Loaded copy/mirror probes = Professional to 4e-10. Tests
  `mesh_copy_timestep1`, `mesh_mirror_timestep1`.

### C46. `step_start()` mesh-operation residuals vs the Professional (F4 T7; studied in F4 T8, 2026-10-10)

- **C46a order — RESOLVED** `e730052`: Professional order on one index
  (pairwise probes, record order irrelevant): macro, copy, move, switch,
  mirror (GNU was switch, move, mirror, copy). Test `mesh_order_timestep1`.
- **C46b `control_mesh_move` — RESOLVED in the GNU; deliberate divergence**
  `c21b62e`: the Professional does not move the nodes; it sets
  `node_start_refined` to the moved current coordinates, once per index.
  Manual 6.217: the mesh is moved. The GNU moves node and
  `node_start_refined`, now once per index (was every step). Owner
  decision 2026-10-10: keep the GNU (DIVERGENCES C46b, class A).
- **C46c `control_mesh_rotate` — PARKED**: the F4 T7 probe was an invalid
  input. Valid probe = the Professional example `test/other/rotate1.dat`
  (2D axisymmetric, swept 180 degrees into 5 prism6): Professional 14
  nodes / 5 prism6; GNU rc 1 ("Length too small of : node"). The GNU has
  no run-time 2D -> 3D conversion, its `mesh_rotate_3d` makes one
  90-degree copy, and it reads `control_mesh_rotate_angle` as a 2D
  rotation instead of the sweep angle. Effort **large**. DIVERGENCES C46c
  (class B). Neither rotate record is guarded once per control yet.
  Owner decision 2026-10-10: BACKLOG item, not scheduled.
- **C46d `mesh_activate1` — RESOLVED** `7aeac4a` + `bc2235e`:
  `mesh_activate_gravity_time` read index 0 only (the input uses 10: full
  gravity from the first step, final vely 0), no bottom-to-top
  interpolation, default method1, a method-2 element got full gravity on
  the reduced stiffness, and the reduced stiffness scaled the assembled
  matrix only (inconsistent Newton). Its target names node 5: the macro
  numbering was transposed (x slowest) vs the Professional (x fastest).
  Now final vely/veliy/sigyy = Pro to 3e-4 and the target passes. Residual
  (partly activated element loaded uniformly): DIVERGENCES C46d (class B).
  Tests `mesh_activate_gravity_time1`, `mesh_macro_numbering1/2`.
- **C46e `control_mesh_switch`** `fb834c6`: two entries are a swap
  (manual 6.244, Professional: `-x -y` swaps; GNU read them as the axis
  order). The Professional leaves `node_start_refined` unswitched; the
  GNU switches it (manual 6.244). Owner decision 2026-10-10: keep the GNU
  (DIVERGENCES C46e, class A). Test `mesh_switch_pair1`.
- **Severity**: **LOW** (open: C46c large, not scheduled; C46b/C46e
  decided 2026-10-10).

### C47. `force11` 3D section moment 0.65 % off the Professional (F4 T9, 2026-10-10) — RESOLVED 2026-10-11 (F4 T12)

- **Corpus**: Professional `other/force11` (circular ring, `post_calcul
  -materi_stress -force`, `solver -matrix_pardiso`).
- **Symptom**: since `-matrix_pardiso` is solved with SuperLU (F4 T9) the
  run finishes in ~10 s, with rc 1: `node_dof_calcul 1 -mom1y_sig`
  5.80223 vs target 5.76474 +- 0.01 (Professional 5.76460; `nory_sig`
  and `shey_sig` targets met: -0.2 shear). Before, the Bi-CG did not
  finish, so the difference was hidden.
- **Root cause**: `mesh.cc` `mesh_extrude()`, quad9 -> hex27 branch: the
  inner planes (layer mid-planes and the layer tops below z1) were
  fresh per-element node blocks, so neighbouring hex27 were connected
  only at z0 and z1 (15600 nodes vs Professional 12000; planes
  2400/3600/3600/3600/2400 vs 2400 each; the ring 7.9x too soft). The
  section moment is equilibrium-based, so it was only 0.65 % off.
- **Fix** (F4 T12, `24301c8`): plane q (1..2nq-1) of source node inod is
  node `next_node + (q-1)*nbase + inod`, created once and shared. After:
  `force11` rc 0, `mom1y_sig` 5.764609 (Professional 5.764612), every
  nodal displacement = Professional to 5.5e-7 (relative to the field
  maximum), section forces within 2e-4. Test
  `extrude_quad9_shared_nodes1` (RED old GNU: rc 1, 84 nodes; GREEN: 75
  nodes, displacements = Professional). Residual: hex27 nodal stresses
  (C57).
- **Severity**: **LOW** (resolved). **Effort**: cheap.

### C57. hex27 nodal stress recovery and `post_point` interpolation differ from the Professional (F4 T12, 2026-10-11) — OPEN

- **Corpus**: Professional `other/force11` (after the C47 fix), test
  `extrude_quad9_shared_nodes1`.
- **Symptom**: with nodal displacements equal to the Professional (5.5e-7)
  the `force11` hex27 nodal stresses differ by up to 3.38 % of the field
  maximum for sigxx/sigyy/sigzz (sigxy 0.9 %). The largest differences
  are at the loaded points (x 0, y +-50.1, every z plane: sigxx -816.80
  vs -789.09, +3.5 %); over the nodes the ratio GNU/Professional spans
  0.987..1.035 with mean 1.0001. Second symptom (same test family):
  a `post_point` placed exactly on a hex27 node returns an interpolated
  displacement that is not the nodal value (`extrude_quad9_shared_nodes1`
  draft, point (2, 1, 0.5): disy -3.51050e-5 vs nodal -3.51802e-5 =
  Professional post_point value).
- **Root cause**: not investigated (Gauss-to-node stress recovery for
  hex27, cf. C42 for quad4; hex27 shape functions / point location in
  `post_point`).
- **Severity**: **LOW**. **Effort**: medium.

### C58. `control_restart`: element history records outside the integration-point dofs not verified (F4 T6c, 2026-10-11) — OPEN

- **Context**: F4 T6c made `control_restart -yes` rebuild the
  integration-point state (`restart.cc` clears `ELEMENT_DOF_INITIALISED`;
  the Professional deletes `element_dof` / `element_intpnt_dof`), tests
  `test-2026/control_restart_reset_ip1/2`.
- **Verified** (probes, Professional 02-08-2026): `element_truss_force`
  is kept by both binaries (1.5 after restart); the undrained pressure
  (`element_intpnt_materi_undrained_pressure`) = Professional (-9.75297).
- **Not verified**: interface IP strain/stress/tension status and
  `element_interface_strain_normal`, `element_contactspring_force`,
  `element_spring_strain`, `element_tendon_strain/_stress`,
  `element_intpnt_materi_stress_effective` (damage). Each is read back as
  history from VERSION_NORMAL and is not reset by `restart()`.
- **`number_timesteps`**: the Professional sets it to 0 on a restart
  (`number_timesteps_total` kept); the GNU has no such record.
- **Next**: one restart probe per family vs the Professional; reset the
  records the Professional deletes.
- **Docs**: `manual-user/` and `manual-developer/control_restart.md`,
  DIVERGENCES C58.

### C59. `group_truss_memory` / `materi_displacement` input errors of the Professional (F4 T6d, 2026-10-11) — RESOLVED 2026-10-11 (F4 T6d)

- **Resolution** (owner decision 2026-10-11): `truss.cc` stops with an
  input error naming both records for an explicit `-updated` /
  `-updated_linear` with `materi_displacement` and for `-total` /
  `-total_linear` without it; `check_error -no` skips both checks (as the
  Professional). sfnet `truss5` adapted (external copy
  `external-downloads/sfnet/extracted/test/truss5.dat` and the omen
  corpus): `-updated` dropped, default `-total_linear`, target = the
  Professional run of the adapted input (node 1 rhside vely -93.8). The
  copy `test/truss5.dat` inside `validation-suite/test.zip` is NOT
  rewritten (archive) and still holds the original. No other corpus input
  hits either check. Tests: `test-2026/truss_memory_updated_displacement_error1`,
  `truss_memory_total_velocity_error1` (expected error).

Original report:


- **Professional 02-08-2026**: `-updated` (or `-updated_linear`) with
  `materi_displacement` stops with "materi_displacement should not be
  initialised"; `-total` / `-total_linear` without it stops with
  "materi_displacement should be initialised" (both suppressible with
  `check_error -no`). The GNU accepts every combination.
- **Owner request** (2026-10-11): make `-updated`/`-updated_linear` +
  `materi_displacement` an input error. Corpus check before adding it:
  sfnet `truss5` (and its copy `test/truss5.dat` in
  `validation-suite/test.zip`) uses exactly that combination (the GNU runs
  it rc 0, the Professional stops). test-2014 `gen1`, `hreltime`, sfnet
  `genera1/2` use `-updated_without_rotation` (another type). No hit in
  test-2026 or the Professional test tree.
- **Next**: owner decision on `truss5`; then the error with a GNU message
  naming both records, plus an expected-error test.
- **Docs**: `manual-user/` and `manual-developer/group_truss_memory.md`,
  DIVERGENCES TRUSS-MEM-ERR.

### C48. Interface tension status written from the wrong buffer (F4 T9, 2026-10-10) — RESOLVED 2026-10-10 (F4 T10)

- **Found by**: the F4 T9 ASan runs on the Professional
  `interface_bar2_hex8`, `interface_bar3_quad8`, `interface_quad8_hex20`
  (stack-buffer-overflow READ in `db()` from `interface_element()`,
  `interface.cc` ~1031), once the `interface_convert()` use-after-free was
  fixed.
- **Root cause**: `element_interface_intpnt_materi_tension_status` is an
  INTEGER record; `interface_element()` fills a `double` buffer allocated
  with `get_new_int()` (released with `free()`, also `rec_stress` /
  `rec_strain`) and PUTs it with `idum[1]` as the integer array and length
  `n_intpnt`: the stored status values are read from the stack (garbage).
- **Fix**: fill a `long int` array and PUT it as integers; release the
  three buffers with `delete[]`. Changes that output record, so it needs
  its own A/B.
- **Severity**: **LOW** (status output only; stresses and strains are
  written correctly). **Effort**: **cheap**.
- **Resolution** (F4 T10, commit `c06fa2d`): `long int` status array
  (`-CLOSED` / `-OPENED`) PUT as integers, buffers released with
  `delete[]`. The record is print-only (nothing reads it back), so only
  that output changes. Values equal to the Professional 02-08-2026 on its
  20 `interface*` inputs (all `-closed`; `interface_bar2_hex8` differs by
  element numbering, see C50) and on the new
  `test-2026/interface_tension_status1` (`-opened` case). Two GNU runs give
  byte-identical `.dbs`.

### C50. `interface_bar2_hex8`: extruded + converted interface mesh differs from the Professional (F4 T10, 2026-10-10)

- **Found by**: the F4 T10 tension-status sweep (Professional `test/other`).
- **Symptom**: after `control_mesh_extrude` (`_n 2 2`) and
  `control_mesh_convert` of a 3D `-bar2` interface the GNU final mesh holds
  element 3 as a `-quad4` (with interface output records of 2 integration
  points) and `-hex8` elements 5-10 with interfaces 9 and 10; the
  Professional has `-hex8` elements 1-6 with interfaces 3 and 6. Both end
  rc 0 (target met).
- **Root cause**: not investigated.
- **Severity**: **LOW**. **Effort**: **medium** (unknown).

### C49. The build uses `-O1` (F4 T9, 2026-10-10) — RESOLVED 2026-10-11 (F4 T13, default `-O2`)

- **Symptom**: `scripts/build_safe.sh` compiles with `-O1`. The F4 T9
  investigation measured about 25 % less run time with `-O2` on the
  solver-bound inputs.
- **To do**: find why `-O1` was chosen (history of `build_safe.sh` and of
  the makefile flags; numerical reproducibility, a past miscompilation, or
  the incremental-build corruption on the NAS) before changing it; then a
  full corpus A/B (`-O2` may change round-off).
- **Study (F4 T11, 2026-10-10; default NOT changed, owner decision)**:
  - **Why `-O1`**: no recorded reason. It comes from the 2026-08-02
    modernisation build command (`VALIDACION-2026-08-02.md`) and was
    copied into `build_safe.sh` when it was created (`73f0e70`, an
    "anti-OOM WSL2" sequential build with a 4 GB `ulimit -v`). The
    upstream makefile Linux targets use `-O2`. No commit or document
    reports a miscompilation or a reproducibility problem at `-O2`;
    `0c2f5e0` (uninitialised `ddum3`) shows UB surfacing with a compiler
    change at `-O1`.
  - **Build** (local, clean, `build_safe.sh` + `TN_EXTRA_FLAGS`, which
    is appended after `-O1`): `-O1` 35 s / largest compiler process 137 MB,
    9 warnings; `-O2` 44 s / 149 MB, 22 warnings; `-O2
    -fno-strict-aliasing` 45 s / 151 MB; `-O3` 49 s / 158 MB, 33 warnings.
    Far below the 4 GB cap.
  - **Speed** (omen, one run at a time, 3 repetitions, medians):
    `large1` 23.4 / 16.9 / 22.2 / 16.2 s and `force11` 14.5 / 8.4 /
    14.0 / 7.7 s (O1 / O2 / O2 -fno-strict-aliasing / O3). Corpus wall
    on the tests finished by both (Pro other, 6 jobs) 731 -> 623 s (O2),
    695 s (O2 nsa), 593 s (O3); sfnet 221 -> 206 / 212 / 190 s. The gain
    comes almost entirely from type-based (strict) aliasing.
  - **Results**: every `.dbs` of the tests finished by both builds is
    byte-identical to `-O1` (O2, O2 nsa, O3; the only difference is the
    2-thread `bicg_csr_threads1`, run-to-run round-off). Peak memory the
    same. Pro other `alternate1`, `cylinder_buckle` finish within the
    120 s timeout at O2/O3.
  - **Blockers**: at O2, O2 nsa and O3 the Pro `interface7`,
    `interface10`, `mesh_cut_2` and the vs `creset_iface` stop with a
    database error (`element_interface_strain_normal`,
    `control_mesh_cut_node_force`), and build_safe gives 177 OK / 4 FALLO
    (`creset_iface`, `bmom1`, `bmom_switch`, `bmom_truss`, `bmom_2d`:
    the cut distance column comes out 0). These are latent stack
    overflows of 1-element scratch buffers given to multi-value `db()` /
    array calls, the C48 class, which the `-O1` stack layout happens to
    tolerate: `print_beam_force_moment.cc` `array_distance(p1, p2, ddum,
    3)` with `ddum[1]`; `print_interface_stress.cc` reads 4 coordinates
    into `ddum[1]` (gcc -Warray-bounds at O2) and the per-intpnt
    `ELEMENT_INTERFACE_FORCE_NORM` / `_STRAIN_NORMAL` into scalars;
    `delete.cc` reads `CONTROL_MESH_CUT_NODE_FORCE` (up to 3 integers)
    into `idum[1]`; `top.cc` reads up to 3 values of
    `control_timestep_iterations_automatic` into an array of 2.
  - **Recommendation**: keep `-O1` for now. First fix this buffer class
    (audit of `db()` / array calls with `idum` / `ddum` scratch buffers,
    then the suite and the corpus under the ASan build), then move to
    `-O2` (not `-O3`: little extra gain, more warnings) with a full A/B.
    Expected gain about 28 % (`large1`) to 42 % (`force11`).
- **Resolution (F4 T13, 2026-10-11, owner decision 2026-10-10)**:
  - **Audit**: `scripts/audit_db_buffers.py` (rerunnable, exit 1 on a
    hit) cross-checks every `db()` read, `get_group_data()`, `array_*()`
    with a constant length and the `geometry()` / `msf_element_rule()`
    output arrays against `data_length` (`database.cc`) and the callee
    contracts. HEAD `dca225a`: 3450 calls, 45 hits at 34 call sites, all
    real or latent overflows, all fixed in `40ae7cd` (no false positive
    left after the script resolved struct members and prototype
    parameters); 38 run-time-record calls reviewed by hand (C54); 273
    unresolved buffers (pointers, array parameters, `&coord[inol*ndim]`).
    Contract and per-site table:
    `manual-developer/db_buffer_contract.md`.
  - **ASan + UBSan** (`-fsanitize=address,undefined -g -O1`): build_safe
    28 failing runs / 18 failed checks -> 1 run (`dofrhside`, C55) /
    1 check; corpus on omen (6 jobs, 360 s, `hard_rss_limit_mb=4000`)
    45 tests with an ASan report at 9 sites -> 2 tests (`dofrhside` C55,
    Pro `ground19_water_under_dam` C56). No UBSan report in either.
  - **`-O1` unchanged**: build_safe 187 OK / 0 FALLO, 226/226 runs; corpus
    `.dbs` identical to HEAD (only the 2-thread `bicg_csr_threads1`
    round-off).
  - **`-O2` A/B** (same sources): every `.dbs` finished by both builds
    byte-identical; `creset_iface`, `bmom*`, Pro `interface7`,
    `interface10`, `mesh_cut_2` run (they failed at `-O2` in F4 T11);
    Pro `alternate1`, `cylinder_buckle`, `ground11_nonsaturated` now
    finish inside the 120 s corpus timeout. build_safe at `-O2` 187 OK /
    0 FALLO, rc 0. Timing in `SEGUIMIENTO-CONVERGENCIA.md` (F4 T13 row).
  - `scripts/build_safe.sh` compiles with `-O2`; the sanitizer targets
    keep appending `-O1`.

### C54. Run-time record names read into a 1-value dummy without a type check (F4 T13, 2026-10-11) — OPEN

- **Symptom**: the `scripts/audit_db_buffers.py --all` REVIEW list. Most
  run-time `db( data_item_name, ... )` reads dispatch on `db_type()` and
  give the matching buffer `DATA_ITEM_SIZE` values, but
  `post.cc` `post_integrate()` / `parallel_post_node()`, `print_el.cc`
  `print_element()`, `repeat.cc` (`control_repeat_save`,
  `control_repeat_until_item`), `data.cc` (arithmetic, spatial resets) and
  `calcul.cc` assume a DOUBLE record and pass `idum[1]` as the integer
  buffer. An input that names an INTEGER record there overflows `idum`.
- **To do**: check `db_type()` (error for an INTEGER record, or read it
  into a long buffer). No corpus input triggers it.
- **Severity**: LOW. **Effort**: cheap.

### C55. `print_unknownsrhside()` reads `node_rhside` with the wrong index (F4 T13, 2026-10-11) — RESOLVED 2026-10-11 (F4 T6a)

- **Symptom**: ASan heap-buffer-overflow READ in `print_un.cc:149`
  (test-2014 `dofrhside`, build_safe and corpus).
- **Cause**: `NODE_RHSIDE` holds `npuknwn` values per node
  (`data_length = npuknwn`) but the loop prints `node_rhside[iuknwn]`
  with `iuknwn = ipuknwn*nder`: the column is wrong for every principal
  unknown after the first when `nder > 1` and the last one reads past
  the record. The fix (`node_rhside[ipuknwn]`) changes the printed
  `*_rhside.N` files, so it needs its own verification (Professional
  `node_rhside` print).
- **Severity**: LOW (print only). **Effort**: cheap.
- **Status**: **RESOLVED 2026-10-11 (F4 T6a)**: `node_rhside[ipuknwn]`.
  Professional 02-08-2026 probe with `derivatives` (one quad4 pushed
  down): Pro `vely_rhside` -1.0416667 (bottom) / +1.0416667 (top),
  `sigyy_rhside` -2.0833333; old GNU 0 / 0 and 4.7e-310 in `disy`; new
  GNU = Pro (= old GNU without `derivatives`). Test
  `dof_rhside_derivatives1` (build_safe file check, RED old / GREEN new).

### C56. Parser tokens longer than `MCHAR`-1 overflow the input buffers (F4 T13, 2026-10-11) — RESOLVED 2026-10-11 (F4 T6a)

- **Symptom**: ASan stack-buffer-overflow (`strlen` in
  `input_skip_comment()`, `input.cc`) on the Professional
  `ground19_water_under_dam`: line 2 has a 141-character URL inside a
  comment; `cin >> str` (`char str[MCHAR]`, 100) writes past `str`
  into `str_tmp`. Results are fine at `-O1` and `-O2` by layout luck.
- **To do**: read tokens with a width limit or into `std::string`
  (11 `cin >> str` sites in `input.cc`) and decide between truncating
  inside comments and an error for an over-long data token.
- **Severity**: LOW. **Effort**: cheap.
- **Status**: **RESOLVED 2026-10-11 (F4 T6a)**: `input.cc`
  `input_word()` reads at most `MCHAR`-1 characters (`setw`) at all 16
  word reads of the parser (11 `cin`, include file, 2 `.run`, 2 gmsh); a
  longer word is split inside a comment (the parenthesis count sees
  every piece) and is an input error (`Word too long`) elsewhere. Tests
  `input_long_comment_word1` (rc 0) and `input_long_word_error1`
  (expected error; old GNU rc 139). Not converted: the two `cin >> str`
  reads of `range.cc` (outside the work unit; manual-developer
  `input_word_length.md`).

### C52. Post-solve material state pass for inertia models (`dynamic3`, F4 T5, 2026-10-10) — RESOLVED 2026-10-11 (F4 T6d)

- **Symptom**: with the ungated pass (TN_C35 prototype) Pro `dynamic3`
  (1D truss, `inertia_apply -yes`, prescribed acceleration 1.234,
  analytic disx 61.7) moved away from the Professional (worst relative
  difference 0.45 -> 0.74, disx).
- **Status**: gated off (`post_solve_state_pass_gated()`: inertia).
  Needed: a pass that updates the inertia terms consistently with the
  final iterate (acceleration/velocity history of the truss), then probes
  + full corpus A/B (widening criterion, C35).
- **F4 T6d attempt (2026-10-11)**: after the C11 fix the `dynamic3` gap
  was the truss memory default / sign, not the pass (c11 = Pro 60.44025,
  identical with the inertia gate removed). Widening the inertia gate
  alone (`top.cc`, materi_acceleration / inertia_apply checks removed) was
  NOT clean: full corpus A/B c59 20b5dc20 -> c52 b4f0f5eb (omen) test-2014
  `iface_gen`, `iface_gen_geom`, `iface_gen_method` exit 0 -> 1
  (interface elements + `inertia_apply -yes`: node 6 disy -0.1009 ->
  -7.417, vely 0.0007 -> 0.551; `iface_gen_quad9` 9 % node_dof, targets
  met); Pro other, sfnet, test-2026 without class change (`taylor3`
  4e-4). build_safe with c52: 223/226 runs. Not committed; the gate stays.
  Next: make the pass consistent for interface elements under inertia
  (or keep interfaces gated) before widening.
- **Root cause of the `iface_gen*` break (F4 T6d, 2026-10-11)**: not the
  pass but the interface Mohr-Coulomb history. `interface.cc` stored the
  ELASTIC trial `kt*gamma_total` in `ELEMENT_INTERFACE_FORCE_TANG`, so the
  friction force had no memory of the plastic slip and did not unload
  after a slip reversal; the one-iteration lag (C35) masked it.
  Quasi-static reversal probe: Professional 02-08-2026 4.000000083,
  analytic 4.0, GNU 5.0. Fixed (owner decision 2026-10-11): the history
  is the clamped force (`b318484`); test-2026
  `interface_mc_slip_reversal1`; test-2014 `iface_mc_slip`/`_1step`/
  `iface_mc_num` original targets restored (CORPUS-STALE-TARGETS §6).
  The inertia part of the gate is removed in a separate work unit only
  if the full corpus A/B is clean.
- **Resolution (F4 T6d, 2026-10-11, `c276931`)**: inertia removed from
  `post_solve_state_pass_gated()` (`top.cc`; slide/contact stay gated,
  C53). RED/GREEN: test-2014 `iface_gen`, `iface_gen_geom`,
  `iface_gen_method` rc 1 with the gate removed and the old interface
  history (node 6 vely -1.70554), rc 0 with both changes. Full corpus A/B
  on omen (interface fix -> gate removed): no class or exit-code change
  except the `alternate1` timeout flake; .dbs differences only in inertia
  models (Pro other `taylor3` 4.1e-4, distance to the Professional
  unchanged: 0.2667 -> 0.2649 / 0.3907 -> 0.3926; test-2014 `iface_gen*`
  up to 0.15, targets met; the rest <= 1.3e-7). build_safe 214 OK /
  0 FALLO, runs 226/226. Remaining interface + inertia gap -> C61.
- **Severity**: LOW. **Effort**: medium.

### C61. Interface Mohr-Coulomb under `inertia_apply`: first-step velocity differs from the Professional (F4 T6d, 2026-10-11) — OPEN

- **Symptom**: probe (two quad4 blocks, quad4 interface kn = kt = 1000,
  phi = 0.785398, c = 0, density 1, `inertia_apply -yes`, volume force
  on the right block, right edge velx -0.1, dt 0.1): after one step
  node 6 vely GNU 0.412250 vs Professional 0.273356; after 10 steps
  -0.105 (GNU, pass on) / -0.098 (GNU, gated) vs -0.0456 (Pro);
  interface shear IP 1 1.65 / 1.84 vs 3.96. Same gap with and without the
  post-solve pass, so it predates C52.
- **Cause**: not established (interface + inertia coupling in the first
  step; candidates: the interface rhs/tangent with the dynamic
  predictor, or the initial acceleration).
- **Severity**: LOW. **Effort**: medium.

### C60. Interface `element_interface_intpnt_stress` normal sign depends on the node order (F4 T6d, 2026-10-11) — OPEN

- **Symptom**: quasi-static interface probe (quad4 interface between two
  blocks, the right block pushed against the left one, kn = 1000,
  normal displacement -0.01): the Professional 02-08-2026 writes the
  normal stress -10 (compression negative) for every node order of the
  interface element (`2 4 9 10`, `9 10 2 4`, `4 2 10 9`); the GNU writes
  +10, +10 and -10. The shear (4.0 / -4.0) and the nodal results agree.
  `element_interface_intpnt_strain` also differs (GNU writes the last
  step values, SEGUIMIENTO `element_interface_intpnt_strain` PENDIENTE).
- **Cause (hypothesis)**: the GNU internal sign convention (compression
  = positive `strain_normal`, normal orientation from the node numbering,
  `group_interface.md` RF-5) leaks into the output record.
- **Severity**: LOW (output record only). **Effort**: cheap/medium.

### C53. Post-solve material state pass for slide models (`slide4`, F4 T5, 2026-10-10) — OPEN, investigation

- **Symptom**: with the ungated pass Pro `slide4` (quad4 on a frictional
  `slide_geometry`, `slide_plasti_friction`, `slide_stiffness`) moved
  away from the Professional (velocity 0.49 relative; before 0.024 sig).
- **Cause (hypothesis)**: the slide forces and friction state are applied
  in `slide()` outside the element loop; the pass re-evaluates the
  material without re-applying them.
- **Status**: gated off (slide and contact data). Needed: call `slide()`
  (and `parallel_contact` for contact) consistently in the pass, then
  probes + full corpus A/B (widening criterion, C35).
- **Severity**: LOW (gated). **Effort**: medium.

### C51. Hardsoil yield read the nodal kappa of the previous evaluation of the same step (F4 T5, 2026-10-10) — RESOLVED 2026-10-10

- **Symptom**: test-2014 `mhardsoil_plast` (one quad4, zero total
  strain, initial sigxx -4) alternated with the number of equilibrium
  iterations: kap 0.0030 / sigxx -1.480 with 1 or 3 iterations, kap
  0.0006 / sigxx -2.47 to -2.49 with 2, 4, 8, 20, 50 (`mhardsoil_gp0_off`
  likewise 0.0015 / 0.00029). The F4 T5 post-solve material state pass
  adds one evaluation and flipped both tests to the odd branch.
- **Root cause**: `plasti.cc` `plasti_rule()` takes `gamma_p` from
  `new_unknowns[kap_indx]`, the interpolated NODAL kappa written by the
  previous material evaluation of the step, instead of the step-start
  kappa plus the increment of its own return mapping.
- **Fix**: `stress.cc` `set_stress()` passes the yield rules a copy of
  the unknowns whose kappa and kappa_shear slots are `old + increment of
  this return` (`pl_unknowns`). Also covers the kappa-hardened
  Mohr-Coulomb variants and the von Mises Nadai hardening (same read).
- **Verification**: no Professional oracle (zero field on the converted
  inputs). Identical results for 1/2/3/4/8 iterations with and without
  the post-solve pass; end state satisfies `sigma = sigma0 - C(E50,nu50)
  eps_p` exactly (HEAD did not), `kap = sqrt(0.5 eps_p:eps_p)`, f ~ 1e-8.
  New tests `test-2026/hardsoil_kappa_iterations1/2`. Retarget (project
  owner 2026-10-10, GNU-calibrated, no Pro oracle): `mhardsoil_plast`
  kap 0.00104991 / sigxx -3.11898 (was 0.0006 / -2.44),
  `mhardsoil_gp0_off` kap 0.000494563 / sigxx -1.56781 (was 0.000293 /
  -1.225).
- **Open, same pattern (not changed)**: `plasti.cc` reads other history
  dofs from `new_unknowns` inside the yield rules: the cap `epp_vol`
  (~221), the cap1 `pc` (~282), the void fraction (~957/987). Each needs
  its own probe (iteration-count sweep) before a change. The nonlocal `f`
  (`fn_indx`) is nonlocal by design.
- **Severity**: MEDIUM (iteration-count dependent results). **Effort**:
  cheap.

---

## D. Summary

| id | corpus | item | severity | effort |
|---|---|---|---|---|
| B1 | sfnet | `examp6` remesh node numbering | MEDIUM | cheap/large |
| B2 | sfnet | `examp14` stale target + runtime | LOW | cheap |
| B3 | sfnet | `examp17` stale target | LOW | medium |
| B4 | sfnet | `examp24` O(N²) mesh mapping | LOW | large |
| B5 | sfnet | `ho_mech2` stale p-refinement input | LOW | medium |
| B6 | sfnet | `ho_mech3` model divergence | LOW | large |
| B7 | sfnet | `ho_mech4` missing strain record | LOW | medium |
| B8 | sfnet | `petsc1` PETSc infrastructure | LOW | large |
| B9 | sfnet | `taylor3` stale input + solver | LOW | medium |
| B10 | sfnet | `condif2/3`, `examp1`, `incnav5`, `strconv1` targets calibrated on the legacy `convection_apply -yes` default (K03) | LOW | cheap |
| B11 | sfnet | `incnav3` target predates the `expansion_volume` inertia density (K04) | LOW | cheap |
| B13 | sfnet | `axisym3/4`, `contac1/2`, `examp12`, `temp4/5`, `visc_pl2` targets calibrated on the old GNU `-follow_material` default (F3) | LOW | cheap |
| B14 | sfnet | `calcul1/2`, `ctvmis2`, `examp18` targets calibrated on the old GNU `-size_tot`/`-size_dev` (F4 T6a) | LOW | cheap |
| B15 | sfnet | `truss2`/`truss3` targets calibrated on the legacy `-updated` truss memory default (F4 T6d; Pro gives the new GNU value) | LOW | cheap |
| C1 | sfnet | `pruebavel` dead/malformed | LOW | cheap |
| C2 | validation-suite | `groundflow_interface` `gvelx` sign | LOW | medium |
| C3 | sfnet | `gforce10` 3D section sign — sign resolved 2026-10-11 (F4 T6a: plot_switch convention, test-2014 input restored, -5000 = Pro); square-face tie kept as deliberate divergence (owner 2026-10-11) | closed | - |
| C5 | — | `contact_target_geometry` index-0 alias — resolved 2026-10-07 (F2 T3a) | MEDIUM | medium |
| C6 | validation-suite | `iface_gen*` generation-coverage gap — resolved 2026-10-11 (F4 T6a: generated elements removed from the inputs) | MEDIUM | cheap |
| C7 | — | Mazars gamma with shear; damage + early plasticity path | MEDIUM | medium |
| C8 | — | default `mesh` `-follow_material` vs Pro `-fixed_in_space` — resolved 2026-10-07 (F3 T2) | MEDIUM | large |
| C9 | — | `-total_linear` + compressibility: Pro time-step dependent (deliberate divergence) | LOW | cheap |
| C10 | Pro `other` | static pore pressure tied to `post_calcul -total_pressure` in the Pro (`excavate1` target) | MEDIUM | cheap |
| C11 | — | `control_materi_dynamic` and implicit truss dynamics vs Pro — resolved 2026-10-11 (F4 T6d) | MEDIUM | medium |
| C12 | — | bicg stop does not abort the GNU on a singular system | LOW | cheap |
| C13 | — | `control_solver_bicg_stop` at the timestep index: Pro reads index 0 only (deliberate divergence) | LOW | cheap |
| C14 | — | `element_normal`: curved quadratic sides, `node_normal`/`element_boundary`/`global_normal` not implemented (zero normal and plane MC direct return resolved) | LOW | cheap |
| C15 | — | plane direct cut-off phi_flow factor M/G vs Pro 2 (deliberate divergence) | LOW | cheap |
| C18 | — | default memory with `materi_displacement`: Pro behaves as `-total_linear` — resolved 2026-10-07 (F3 T2) | MEDIUM | large |
| C19 | — | explicit `-updated_linear` on `-follow_material`: GNU reference gradient vs Pro current — resolved 2026-10-07 (F3 T2) | LOW | medium |
| C20 | — | `-updated_jaumann` not implemented — resolved 2026-10-08 (F3 T2) | LOW | medium |
| C21 | — | `control_materi_updated_apply -yes`: Pro ignores its manual (deliberate divergence) | LOW | cheap |
| C22 | Pro `other` | `hypo12` marginal sigyy target miss after the default-memory change | LOW | cheap |
| C25 | — | uninitialized-buffer audit follow-ups (`pef_vel_prev` realloc, NODE_EIGEN print, split/beam/interface/mpc/nonloc side findings) | LOW | cheap |
| C26 | — | `control_print_dof -separate_index`: Pro writes one file per dof, GNU one `dof.N` | LOW | cheap/medium |
| C27 | — | FRD/gmsh layout vs Pro (PSTEP per block, `ALL` component, names, extra blocks/fields, multi-component gmsh fields) | LOW | medium |
| C29 | Pro `other` | `group_condif_conductivity` per-direction values (condif11) not parsed | LOW | cheap/medium |
| C30 | Pro `other` | `mohr_coul_direct4`: Pro does not enforce the interface MC strength (deliberate divergence, GNU = manual 6.632) | LOW | — |
| C31 | Pro `other` | `tet10` nodal quadrature only degree 1 (`tet10_2` gravity) — resolved 2026-10-09 (F4 T4: Keast 5-point + LSQ recovery + HRZ lumping) | MEDIUM | medium |
| C32 | Pro `other` | `axisym9` (stale row, already = Pro), `distri2` (random generator): not GNU defects — closed 2026-10-09 (F4 T4) | LOW | — |
| C33 | — | `_visco` direct cut-offs: `DTIME` from the previous step, no visco on the principal tension cap — resolved 2026-10-09 (F4 T3) | MEDIUM | cheap/medium |
| C34 | — | direct cut shear plastic strain: Pro twice the tensor value (deliberate divergence) | LOW | — |
| C35 | — | one-iteration lag of the whole material state (K02, `spring4` 1 iteration) — resolved 2026-10-10 (F4 T5, post-solve material state pass); inertia/slide/contact gated (C52, C53) | MEDIUM | medium/large |
| C36 | validation-suite | `young_power` + `young`: Pro adds both stiffnesses, manual says law = modulus (K06) — decided 2026-10-06: deliberate divergence (GNU = manual), no code change | LOW | cheap |
| C37 | Pro `other` | planar GNU beam vs Pro 3D frame beam (`trubea2/3`, K17) — parked, design | MEDIUM | large |
| C38 | Pro `other` | `control_support_edge_normal_stiffness_freeze` parse-only (support_edge_normal_3, 1e-9 effect) | LOW | cheap/medium |
| C39 | Pro `other` | K17 residuals: `truss8` birth node resolved 2026-10-06 (death velocity reset open), `update2` -updated robustness, `truss13` truss_distribute_mpc, zero-force truss first iteration | LOW | cheap/large |
| C40 | Pro `other` | element-target contact residuals after F2 T3b: no Pro check "target element has all nodes on the boundary", hex8 target sides never detected, node_boundary generated only for element-target contact | LOW | cheap/medium |
| C41 | Pro `other` | `mohr_coul_direct7`: `-total_linear` + velocity `mpc_node_number` tie — tie resolved 2026-10-09 (F4 T2); residual disy/sigzz difference | LOW | medium |
| C42 | Pro `other` | quad4 nodal stress recovery: natural bilinear (GNU) vs physical linear least squares (Pro); distorted/axisymmetric quad4 (`axisym2` sigzz 10760.09 vs 10759.22) | LOW | medium |
| C43 | Pro `other` | memory hogs of the Bi-CG solver (`force11`, `large1`, `ground19`, `mesh_activate1`, `foundation_consolidation_settlement`): resolved 2026-10-09 (F4 T7 A dense n x n matrix, B lazy band retry, C macro per step); solver speed resolved 2026-10-10 (F4 T9: `-matrix_pardiso` -> SuperLU, CSR operator, recovery table); `foundation...` still Bi-CG bound | MEDIUM | cheap/medium |
| C44 | Pro `other` | one-shot mesh records (`control_mesh_macro`/`_copy`/`_mirror`/`_switch`) re-applied on every time step — resolved 2026-10-09 (F4 T7 C) | MEDIUM | cheap |
| C45 | — | `control_mesh_copy`/`_mirror` after the first time step: `node_dof_start_refined` error — resolved 2026-10-10 (F4 T8 `0ca144c`) | LOW | cheap |
| C46 | Pro `other` | `step_start()` mesh residuals — F4 T8: order, move once, switch pair, `mesh_activate1` resolved; rotate sweep parked, not scheduled (C46c, large); move/switch `node_start_refined` deliberate divergences (A, 2026-10-10) | LOW | large (C46c) |
| C47 | Pro `other` | `force11` 3D section moment `mom1y_sig` 5.80223 vs Pro 5.76460 — resolved 2026-10-11 (F4 T12: quad9 -> hex27 extrusion shares the inner-plane nodes) | LOW | cheap |
| C49 | — | build at `-O1` — resolved 2026-10-11 (F4 T13: 34 scratch-buffer call sites fixed, default `-O2`, results identical) | LOW | medium |
| C50 | Pro other | `interface_bar2_hex8` extruded/converted interface mesh numbering | LOW | medium |
| C51 | validation-suite | hardsoil `gamma_p` read the nodal kappa of the previous evaluation (iteration-count oscillation, `mhardsoil_plast`) — resolved 2026-10-10 (F4 T5); cap/cap1/void reads open | MEDIUM | cheap |
| C52 | Pro `other` | post-solve state pass for inertia models (`dynamic3`) — resolved 2026-10-11 (F4 T6d: interface history fix + inertia removed from the gate) | LOW | medium |
| C53 | Pro `other` | post-solve state pass for slide/contact models (`slide4`) — gated off, investigation | LOW | medium |
| C54 | — | run-time record names read without a type check into a 1-value dummy (`post.cc`, `print_el.cc`, `repeat.cc`, `data.cc`, `calcul.cc`) | LOW | cheap |
| C55 | validation-suite | `print_unknownsrhside()` indexes `node_rhside` (npuknwn values) with `iuknwn` (heap read past the end, wrong column when nder > 1)  — resolved 2026-10-11 (F4 T6a)| LOW | cheap |
| C56 | Pro `other` | input tokens longer than `MCHAR`-1 (`cin >> str`) overflow the parser buffers (`ground19_water_under_dam`, 141-char URL in a comment)  — resolved 2026-10-11 (F4 T6a)| LOW | cheap |
| C57 | Pro `other` | hex27 nodal stress recovery (`force11` up to 3.5 % at the loaded points, displacements equal) and `post_point` interpolation on a hex27 node | LOW | medium |
| C58 | Pro `validation` | `control_restart`: interface/spring/contact-spring/tendon/damage element history records not verified after the F4 T6c integration-point reset (`number_timesteps` record absent) | LOW | cheap/medium |
| C59 | sfnet | `group_truss_memory -updated` + `materi_displacement` (and `-total*` without it): input error in the Pro — resolved 2026-10-11 (F4 T6d) | LOW | cheap |
| C61 | — | interface Mohr-Coulomb + `inertia_apply`: first-step velocity differs from the Pro (0.412 vs 0.273), independent of the pass | LOW | medium |
| C60 | — | interface `element_interface_intpnt_stress` normal sign depends on the node order (Pro: compression negative always) | LOW | cheap/medium |

Forty-one items remain open: four MEDIUM and thirty-seven LOW (C23 and C24 resolved 2026-10-05; C25, C26, C27 and B10 added; F1 T3 / K04 2026-10-06 added B11, C28-C32; F1 T6 resolved C28, added and resolved B12; F2 T1 added C36; F2 T2 added C37, C38, C39; F2 T3 resolved C5 and added C40; F3 T2 resolved C8, C18, C19, C20 and added B13, C41; F4 T2 fixed the C41 tie and left a residual; F4 T3 resolved C33; F4 T4 closed C32, added C42 and resolved C31; F4 T7-T8 added C43-C46; F4 T9 resolved the C43 solver speed, decided C46b/C46e and added C47, C48, C49; F4 T10 resolved C48 and added C50; F4 T5 added C52, C53 and resolved C51; F4 T13 resolved C49 and added C54-C56; F4 T12 resolved C47 and added C57; F4 T6a resolved C6, C55, C56 and the C3 sign, the C3 square-face tie becoming a deliberate divergence (owner 2026-10-11), and added B14; F4 T6c fixed the `control_restart` integration-point reset and added C58; F4 T6d resolved C11 and C59, fixed the interface history and resolved C52, and added B15, C60 and C61). No open issue is CRITICAL or
HIGH. The validation-suite corpus has no open failure — its 233 real tests
pass (`materi_direct_mc` input corrected, C17 resolved). Every remaining failure is a stale or malformed test, an out-of-scope model
or infrastructure limit, a narrow post-processing divergence, a corpus-coverage
gap, or a performance cost. The source defects found during the cycle were fixed
and validated against the corpora. The remaining post-convergence improvements
(items b–f7) are proposals in `IMPROVEMENT-ANALYSIS.md`; item (a), the post-solve
stress pass hoist, is the only one already implemented.
