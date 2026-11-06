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

### B15. sfnet `truss2`/`truss3` targets calibrated on the legacy GNU truss memory default (F4 T6d, 2026-10-11) — RESOLVED 2026-10-11

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
- **Status**: **RESOLVED 2026-10-11 (`94b0254`)**: both local and omen
  sfnet copies now declare `group_truss_memory 0 -updated` (`truss3` also
  uses `materi_velocity_integrated`), as the Professional copies; GNU and
  the Professional give the original target (0.1, 0.1), exit 0.
- **Severity**: LOW. **Effort**: cheap.

### B16. sfnet `examp3` target calibrated on the nodal-kappa defect (F4 T6f WU2, 2026-10-12) — RESOLVED 2026-10-12

- **Corpus**: sfnet 197 (external). `examp3` (2014 ancestor of the
  Professional `validation_3`) targets `post_node_result 1 -vely` -124.8
  ±1, a 2014 GNU value that encoded the nodal-kappa defect (hardening
  read the h-lumped nodal kappa back at the points). With the
  integration-point kappa of F4 T6f WU2 GNU gives -126.431 (exit 0 -> 1);
  the Professional 02-08-2026 on the translated input gives -126.374.
- **Status**: **RESOLVED 2026-10-12** (owner decision): local and omen
  sfnet copies recalibrated to -126.374 ±1 with a dated header comment
  (CORPUS-STALE-TARGETS §8); GNU exit 0.
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

### C3. `gforce10` — 3D section face/arm sign limitation — RESOLVED 2026-10-11 (F4 T6a: sign; square-face tie class A C3)

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

### C5. `contact_target_geometry` alias only reads index 0 — RESOLVED 2026-10-07 (F2 T3a)

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

### C7. Mazars damage: gamma with shear, and damage + early plasticity — RESOLVED 2026-10-12 (deliberate divergence, class A)

- **Corpus**: none (own hex8 / bar probes, `materi_damage_effective_stress.md`).
- **Symptom**: (a) with a compressive `sigxx` and a same-sign `sigxy`/`sigxz`
  the Mazars `dam` differs from the Professional 02-08-2026 (e.g. 0.10809 vs
  0.18752; with the opposite shear sign 0.26351 vs 0.39352); (b) a bar with
  `group_materi_plasti_vonmises` yielding before the damage threshold gives a
  different plastic strain path (`eppxx` 1.879e-4 vs 1.368e-4).
- **Root cause (F4 T6f, 2026-10-12)**: not gamma. The Professional uses a
  non-standard equivalent strain: it accumulates the positive principal
  strains in ascending order as s = 2 s + eps^2 (weights 4, 2, 1 instead
  of 1, 1, 1), eps_eq = sqrt(s) times gamma; this form reproduces 14
  Professional probes to 7e-14. The GNU follows the manual 2.2.5,
  eps_eq = sqrt(sum <eps_i>+^2). (b) is the same defect (the plastic path
  follows from the different damage). Both codes use the same gamma.
- **Resolution**: owner decision 2026-10-12: keep the GNU (manual), class A
  row C7 in `DIVERGENCES-PROFESSIONAL.md`, SEGUIMIENTO "Descarte de
  features". Guard tests (`ae0691b`): test-2026
  `damage_mazars_eq_two_positive1` (hex8 diag(2e-4,-3e-4,1e-4), dam
  0.28610207515; Pro 0.36409) and `damage_mazars_vonmises1` (dam
  0.925372225695, eppxx 1.87863e-4; Pro 0.96946 / 1.368e-4).
- **Note (no fix planned)**: the triaxiality factor gamma is not objective
  in either code: `damage_mazars()` (`damage.cc:84-95`) builds it from the
  xx, xy, xz row of the stress matrix instead of the principal stresses,
  so two reflected states (xx -5e-4, xy -/+5e-4, zz 3e-4) give different
  `dam`. It only matters with a compressive xx combined with xy or xz
  shear. Shared with the Professional, so not a divergence; recorded as a
  latent model issue (owner decision 2026-10-12: note only).

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

### C12. bicg "stop" does not abort the GNU run on a singular system — CLOSED 2026-10-17 (superseded, F4 T6f LOW queue)

- **Closure (2026-10-17)**: re-measured on the same probe (elasti7 mesh
  without supports, `solver -matrix_iterative_bicg`): the Professional
  02-08-2026 does NOT stop (rc 0, `.dbs` with displacements ~5e9, no
  message); the GNU since C77 prints `Warning: iterative solver
  breakdown, system solved with the direct SuperLU` and also finishes
  with rc 0 (displacements ~1e10). Both codes leave a singular model
  undetected; no divergence remains, the C77 guarded LU fallback
  replaced the per-iteration breakdown error. Probe
  `scratchpad/lowq/c12/`.

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
- **WU1b (F4 T6f K21, 2026-10-12)**: the principal tension cap is lifted
  too (`group_materi_plasti_tension_direct 0 1.e10`), target `sigxy`
  1 +- 1e-6 (GNU exactly 1 every step); same rewrite for `cd_method` and
  `cd_geom` (all `vely` nodes). Professional 02-08-2026: exactly 1 with
  `-updated_linear` (also the old input); with `-updated` it gives
  2.95 .. 18.69 at t = 1 and so does the GNU with `-updated` (18.6915
  vs Pro 18.6915; `cd_method` 18.0374 vs 18.0375). The "Pro 18.0 vs 1.0"
  of the WU1b investigation was a translation artifact
  (`-updated_without_rotation` has no Professional equivalent; on these
  shears it behaves as the Pro `-updated_linear`, both sigxy 1, not as
  `-updated`), NOT an MC-direct divergence.

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

### C22. `hypo12` (Pro `other`) marginal target miss after R10 — DEFERRED (class B, F4 T6f LOW queue, owner plan 2026-10-17)

- **Deferred (2026-10-17)**: a 1.02e-3 miss for a 1e-3 tolerance on a
  Professional file; revisit with the Masin kernel calibration (K14).

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

### C25. Uninitialized-buffer audit follow-ups (F0 T6) — item 1 RESOLVED 2026-10-17 (`efff819`); rest DEFERRED (class B, F4 T6f LOW queue, owner plan 2026-10-17)

- **Status (2026-10-17)**: item 1 fixed (`efff819`: `pef_vel_prev`
  zeroed after every reallocation). The `print_beam_force_moment.cc`
  work array (now `work[3]`) and the `interface.cc` `free()` are already
  gone. Deferred: item 2 (`print_te.cc`/`print_pl.cc` NODE_EIGEN reuse,
  probably unreachable), item 3 (`split.cc` geometry output in the zero
  fallback) and the `mpc.cc` / `nonloc.cc` side findings: outside the
  LOW-queue edit surface, none observed in the corpora.

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

### C35. One-equilibrium-iteration lag of the whole material state (audit K02) — RESOLVED 2026-10-10 (F4 T5), gated families resolved 2026-10-11 (C52, C53)

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

### C37. 3D beam: the GNU beam is planar, the Professional beam is a full 3D frame element (audit K17 `trubea2/3`) — RESOLVED 2026-11-02

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
- **Status**: **RESOLVED 2026-11-02 (F5 T2 V4d)**: the 12-dof local-axes
  beam `beam_local()` of F5 T2 V4c is now used for every 3D beam, with the
  Professional default local z (owner decision "Me gustaría tener vigas 3D
  nativas"); the planar 3D path was removed. `trubea2` / `trubea3` = the
  Professional (node_dof and node_rhside, disz 0.0333, roty -0.005 /
  rotx 0.005); the sfnet `trubea2` target (`node_rhside 1 -roty` +1e-2)
  carries the old planar sign, the Professional gives -1e-2 too.
  Test `beam_3d_default_axes1` (14 directions = Professional probe `obl`).

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
  carries no force in the GNU (Pro 5.0). Since C77 (2026-10-15) the
  resulting singular system is no longer accepted on the old absolute
  floor: the run stops with an explicit breakdown error after SuperLU
  and the band LU both refuse the matrix (rc 1 before and after).
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
  targets are never detected (code reading, not measured on the Pro).
- **(3) Resolved 2026-11-02 (F5 T2 V3)**: `mesh_boundary -yes` (6.816) now
  fills `node_boundary` and `element_boundary` (`area.cc`
  `mesh_boundary_set()`, renamed from `contact_node_boundary_set()`). The
  Professional fills them in every run and marks every element: the GNU
  keeps the manual (default -no, only elements with a free side), owner
  decision 2026-11-02, DIVERGENCES MESH-BOUNDARY-DEFAULT (class A) and
  ELEMENT-BOUNDARY-ALL (class D).
- **Severity**: **LOW**. **Effort**: **cheap** (1, 2).
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
### C42. quad4 nodal stress recovery: natural-coordinate bilinear interpolation (GNU) vs physical-coordinate linear least squares (Professional) (F4 T4, 2026-10-09) — RESOLVED 2026-10-12 (F4 T6f K21 WU1)

- **Resolution**: physical-coordinate linear fit per element for the 2x2
  quad4/quad9 (`miscel.cc` `recovery_fit_factor()`); `axisym_quad4_ring1`
  sigzz 10759.2155 = Professional. The `cd_method` flip of the F4 T4
  prototype was first avoided by gating the fit to non-convected stress
  dofs (`stress_recovery_lagrange`); WU1b (2026-10-12) removed the gate:
  the old `cd_method` was ill-posed (`-ra 1 4 -ra` is a node list, rigid
  rotation, target met by chance) and was rewritten as a well-posed
  shear; the convected stress path is C72.

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

### C57. hex27 nodal stress recovery and `post_point` interpolation differ from the Professional (F4 T12, 2026-10-11) — RESOLVED 2026-10-13 (recovery by C70; `post_point` by F4 T6f C57; curved-element location -> C78)

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
- **Update 2026-10-13 (F4 T6f C70)**: the nodal stress part is RESOLVED:
  the Professional fits the hex27 points in physical coordinates (C70);
  `force11` nodal sigxx/sigyy/sigzz/sigxy now = Professional to 7e-7 of
  the field maximum (before 3.38 %). Open: the `post_point` part.
- **Investigation 2026-10-13 (F4 T6f C57)**: the
  `post_point` part is a GNU extrusion defect, not an interpolation
  one. The quad9 -> hex27 branch of `mesh_extrude()` (`mesh.cc`, the
  shared-plane node loop added by C47 in `24301c8`) writes `NODE` with
  the plane z but copies every other node-class item, including
  `NODE_START_REFINED`, from the source node: the 45 mid / inner plane
  nodes of the witness keep z = 0 in `node_start_refined` (the linear
  branch sets it explicitly). `point_el()` reads `NODE_START_REFINED`, so
  it locates the point in a collapsed element: GNU (2, 1, 0.5) disy
  -3.51050e-5, (1, 0.5, 0.5) -1.13556e-5, (0.5, 0.5, 0.5) -3.25367e-6;
  the exact hex27 interpolation of the GNU nodal dofs gives -3.518023e-5,
  -1.165094e-5, -3.574477e-6 = Professional post_point to 1e-14.
- **Fix** (F4 T6f C57, `cf37de3`): the quad9 branch rewrites
  `NODE_START_REFINED` of every new plane node with `coords[2] = zz`, as
  the linear branch does. The other node-creating mesh branches (mirror,
  copy, rotate 2d/3d, move, switch, convert_quad8/hex20, linear extrude)
  already transform it (checked). Test `test-2026/post_point_extrude_hex27_1`
  (four post points, Professional targets, absolute 1e-11): RED on HEAD
  `8b36a34` rc 1 (disy -3.51050e-5 / -1.13556e-5 / -3.25367e-6 /
  -1.78447e-5), GREEN rc 0 (-3.5180218e-5 / -1.1650940e-5 / -3.5744774e-6
  / -1.8052922e-5; the first point is 9.6e-12 off because the
  `point_el` descent stops about 1e-8 iso short).
  Corpus A/B omen `c57_runs` (958 runs): no classification change; 878
  of 882 non-timeout `.dbs` identical (masking `time_calculation` /
  `time_at_start`); `force11` and `extrude_quad9_shared_nodes1` differ
  only in `node_start_refined` (7200 / 45 inner-plane nodes, now = `node`;
  every result record identical), the new test, and `bicg_csr_threads1`
  (2-thread round-off).
- **Second finding (curved elements, -> C78, deliberate divergence)**: on a curved
  quad9 (2D) and a distorted hex27 the GNU `post_point` is the exact
  isoparametric interpolation (2e-12 / 2e-8 of the nodal dofs; the
  descent stops ~1e-8 short), while the Professional 02-08-2026 locates
  the point with the CORNER nodes only (bilinear / trilinear inversion,
  iso within ~3e-6) and then interpolates with the full quadratic shape
  functions: a `post_point` on the quad9 centre node returns velx
  9.437e-7 for a nodal 1.898e-7. The GNU is analytically right there.
- **Severity**: **LOW** (resolved). **Effort**: cheap.

### C58. `control_restart`: element history records outside the integration-point dofs not verified (F4 T6c, 2026-10-11) — interface/spring VERIFIED 2026-10-17; rest DEFERRED (class B, F4 T6f LOW queue, owner plan 2026-10-17)

- **Status (2026-10-17, `ef9af63`)**: restart probes vs the Professional
  02-08-2026: the interface history (normal/shear -1.5 / 1.5 after the
  restart, -0.5 / 0.5 in a fresh run) and the spring history (0.15 vs
  0.05, `-total_linear` and `-updated`) are KEPT by both programs, so
  nothing is reset; pinned by `test-2026/control_restart_interface_history1`
  and `_spring_history1`. Deferred: contact-spring, tendon and damage
  effective-stress probes, and the `number_timesteps` record (output
  only: the Professional sets it to 0 on a restart, the GNU has no such
  record).

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

### C54. Run-time record names read into a 1-value dummy without a type check (F4 T13, 2026-10-11) — RESOLVED 2026-10-17 (`8edaf3b`)

- **Resolution (F4 T6f LOW queue Q4)**: `post.cc` (`post_node`,
  `post_integrate`), `print_el.cc` and `repeat.cc` (`control_repeat_save`,
  `control_repeat_until_item`) check `db_type()` and stop with `Error: a
  DOUBLE PRECISION data item is required here, not the INTEGER item: <name>`.
  The `data.cc` and `calcul.cc` sites were already safe (type checks or
  internal DOUBLE records). Expected-error tests
  `test-2026/post_node_integer_item_error1`, `post_integrate_*`,
  `control_repeat_save_*`, `control_repeat_until_*` (RED rc 0 / endless
  repeat). DIVERGENCES C54-TYPE (class C: the Professional returns 0 for
  `post_node` / `post_integrate`).

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
  `input_word_length.md`). **Closed 2026-10-12 (F4 T6f ranges)**:
  `range_scan()` reads every word through `input_read_string()` (bounded
  `input_word()`, defines, include files); no raw `cin` read is left in
  `range.cc`. Cosmetic note: the Professional manual calls the list
  `-ra 5 4 8 -ra` (mesh_activate_gravity_element) "a number range"; it
  is a list (manual pages `ranges.md`).

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

### C76. `prism6` `group_integration_points -minimal` ignored (F4 T6f C75, 2026-10-13) — RESOLVED 2026-10-16 (F4 T6f C67 remainder + C76, `dc4a59a`)

- **Symptom**: with `group_integration_points -minimal` the Professional
  02-08-2026 integrates a prism6 with 2 points, the triangle centroid
  (1/3,1/3) at the two Gauss levels zeta = 1/2 +- 1/(2 sqrt 3)
  (`element_intpnt_npoint` 2, `element_intpnt_iso_coord`; manual 6.622:
  "-minimal ... number of nodes in the direction minus 1"). The GNU
  `pol()` (`polynom.cc`) reads the record but always uses its 6-point
  rule (C75).
- **Next step**: add the 2-point branch to `pol()`; measure how the
  Professional recovers 2 points to the 6 nodes (`recovery_fit_kind` in
  `miscel.cc` only accepts npoint 6 for prism6; a linear `1,x,y,z` fit
  is not defined on 2 points), then a dedicated test.
- **Severity**: LOW (non-default option). **Effort**: cheap/medium.
- **Resolved** (2026-10-16, `dc4a59a`): `pol()` honours `-minimal` with
  the Professional rule (centroid x the 2 Gauss levels, weight 1/4 each;
  a 20x2x2 elastic cantilever's Professional nodal displacements equal a
  2-point FE model to 1.5e-13). Recovery: per element
  `sum_g h_n(g) sig_g / sum_g h_n(g)` (Professional residual 5e-13; the
  linear fit is 6e-2 off), arithmetic mean (recovery kind 4,
  `miscel.cc`). Analytic check of the rule: 7 zero eigenvalues per
  element (1 spurious mode) but no global zero mode in any cantilever
  mesh tried; patch test exact; tip-load bending error -18.0 / -7.1 /
  -2.8 % (10x1x1 / 20x2x2 / 40x2x4) against -34.7 / -13.1 / -4.85 % of
  the 6-point rule, so the Professional rule is sound for an explicit
  under-integration option. The default stays the 6-point rule (C75;
  `pol()` now sets the prism6 default to `-maximal` before reading the
  record). Tests `test-2026/recovery_prism6_minimal1`,
  `prism6_minimal_cantilever1`.

### C75. `prism6` integration points differ from the Professional (F4 T6f C70, 2026-10-13) — DECIDED 2026-10-13 (deliberate divergence, class A)

- **Symptom**: the Professional 02-08-2026 `element_intpnt_coord` of a
  prism6 are the triangle CORNERS at the two Gauss levels zeta = 1/2 +-
  1/(2 sqrt 3) (element (0,0,0)-(1,0,1): points (0,0,0.2113),
  (1,0,0.2113), (0,1,0.2113), ...); the GNU `pol()` (`polynom.cc`, rule of
  2026-09-07) uses the interior degree-2 triangle points (2/3,1/6,1/6) at
  the same levels (only the zeta levels were checked then). On a
  distorted prism with a non-linear field the nodal sigzz differs by
  1.6 % of the scale (1.9247 vs 1.8713); with a field in the prism6 span
  `1,x,y,z,xz,yz` both are exact (`test-2026/recovery_prism6_linear1`).
- **Measured** (F4 T6f C75, 2026-10-13): the corner points x 2 Gauss
  levels are used by default, with `-normal`, `-maximal` and with
  `inertia_apply -yes`; `-minimal` uses 2 centroid points (C76). Weights
  1/12 each: the nodal displacements of a 20x2x2 prism6 cantilever (E 1,
  poisson 0.3, tip deflection prescribed) equal an independent corner-rule
  FE model to every printed digit (midspan node 11 disz -0.3121926845;
  GNU = interior rule, -0.3123428787). The element volume does not tell
  the rules apart (detJ is linear in the triangle). Manual 6.622 is silent
  on the triangle points of a prism. prism15 was not re-probed (its
  21-point rule was already measured equal to the Professional).
- **Accuracy**: both rules give a full-rank stiffness (6 rigid modes) and
  pass the patch test (internal forces equal to 1e-16 on a distorted
  prism). Cantilever 10x1x1, end load, triangles in xy, against Timoshenko
  -4031.2: corner (Professional) -56.0 / -26.0 / -9.95 %, interior (GNU)
  -34.7 / -13.1 / -4.85 % for meshes 10x1x1 / 20x2x2 / 40x2x4; triangles
  in xz: -80.0 / -51.2 / -21.8 % against -79.6 / -50.2 / -21.1 %. Prescribed
  tip deflection 20x2x2: reaction 3.357e-4 (corner) / 2.861e-4 (interior),
  analytic 2.481e-4. The vertex rule integrates L_i^2 as A/3 instead of
  A/6, which stiffens the element.
- **Decision** (project owner, 2026-10-13): keep the GNU interior rule as
  a deliberate divergence (DIVERGENCES class A row C75, SEGUIMIENTO
  Descarte). No code change. Probes: scratchpad `c75` of the session
  (`m/`, `beam/`, `rank.py`, `beam.py`, `beam2.py`, `gbeam.py`).
- **Severity**: closed (deliberate). **Effort**: -.

### C74. `validation_1` prints "Bi-CG iteration broke down" about 600 times (F4 T6f C72, 2026-10-13) — RESOLVED 2026-10-13 (F4 T6f C74, message part; solver part -> C77)

- **Symptom**: the GNU `validation_1` channel flow (`validation1_steady1`,
  quad4 with p- and h-refinement, `convection_apply -yes`, `inertia_apply
  -no`, t = 0..3, dt 0.01) printed "Error: the Bi-CG iteration broke down."
  plus "Probably you have a wrong input file" 599 times (599 of its 600
  solves). The run ended rc 0 with the correct steady velocity.
- **Cause**: `so_bicg.cc` Bi-CG breakdown test `|r1.r2| < 1e-16` or
  `|d2^T A d1| < 1e-16`, ABSOLUTE. Near the steady state `|r0|^2` is
  1e-12 .. 1e-8, and the test fires on healthy recurrences
  (`|r1.r2|/(|r1||r2|)` 1e-3 .. 1e-4, `|dAd|/(|d2||Ad1|)` 1e-2 .. 1e-4;
  493 events dAd, 66 r1.r2, 40 both). `so.cc` then re-solves with the
  direct band LU (`solve_band_retry`), exact: HEAD = SuperLU
  (-0.0571954 at t = 3). Not limited to validation_1: 38 of 956 corpus
  runs (Professional `test/other` 330, sfnet 197, test-2026 185,
  test-2014 244) printed it; largest ground11_phreatic_level 2983,
  ground11_nonsaturated 2492, sfnet ho_mech1 613, validation1_steady1
  599, stokes_group1/2 199, vajont_slope_fs 150.
- **Professional 02-08-2026**: default solver PARDISO (direct); with
  `solver -matrix_iterative_bicg` it prints no breakdown on validation_1
  either (same velocity).
- **Resolution** (owner decision 2026-10-13, option A): a breakdown that
  the LU retry solves prints ONE line per solve, "Warning: Bi-CG
  breakdown, system solved with the direct band LU"; the full error block
  (`solve_iterative_bicg_breakdown_report()`) is printed only when the LU
  retry fails too, or on the displacement-primary second solve (no
  retry); `solver_bicg_stop -no` keeps its partial-iterate warning. No
  numerical change: corpus A/B on omen (HEAD `46104a3` vs the fix, 957 runs, `c74b_runs`) has no classification change and 881 of the 882 non-timeout `.dbs` are identical after masking `time_calculation` / `time_at_start` (the other, `bicg_csr_threads1`, runs on 2 threads and differs between two HEAD runs at round-off); the runs with the error block drop from 38 (+ the new test) to 0 (`tslv_bsym` keeps its `solver_bicg_stop -no` warning). Test test-2026
  `bicg_breakdown_lu_retry1` (build_safe: rc 0, the warning, no "wrong
  input file"; old GNU: 9 error blocks).
- **Not done** (option B, -> C77): a scale-relative breakdown test.

### C78. `post_point` location in curved quad9 / hex27: the Professional uses the corner nodes only (F4 T6f C57, 2026-10-13) — DECIDED 2026-10-13 (deliberate divergence, class A)

- **Probes** (2026-10-13, Professional 02-08-2026): two curved quad9
  (2D, midside and centre nodes off the straight positions, 7 post
  points) and one randomly distorted hex27 (5 points); the nodal dofs of
  both binaries agree to 1e-12.
- **GNU**: exact isoparametric location (full quadratic map): the
  `post_point_dof` equals the exact interpolation of the nodal dofs to
  2e-12 (quad9) / 2e-8 (hex27, the descent stops about 1e-8 short); a
  post point on a node returns the nodal value.
- **Professional**: locates the point with the CORNER nodes only
  (bilinear / trilinear inversion) and interpolates with the full
  quadratic shape functions: its values match that model to about
  1e-6 relative; the iso coordinates backed out of its values sit
  within about 3e-6 of the exact corner inversion (its own iteration
  tolerance). A post point on the curved quad9 centre node: velx
  9.437e-7 vs nodal 1.898e-7; hex27 interior point disx -6.5709e-3
  vs -5.9265e-3 (11 %).
  Same family as TET10-PPOINT (tet10: Professional linear from the
  corners, deliberate divergence).
- **Decision** (project owner 2026-10-13): deliberate divergence, class
  A. The GNU keeps the exact isoparametric inversion (a post point on a
  node returns the nodal value); the Professional corner-node location
  is not adopted. Straight-edged elements with centred mid nodes agree
  in both codes.
- **Optional (LOW)**: the `point_el` finite-difference descent stops
  about 1e-8 iso units short; an analytic-Jacobian Newton polish (as in
  `mpc.cc` `mpc_weights_polish()`) would make post points exact to
  machine precision (prism15 1.000000261 vs 1.0, the C57 test point 1
  9.6e-12).
- **Severity**: LOW (decided). **Effort**: -.

### C77. Bi-CG scale-relative breakdown test and relative convergence floor (F4 T6f C74, 2026-10-13) — RESOLVED 2026-10-15 (F4 T6f C77; non-symmetric algorithm choice pending)

- **Symptom**: the Bi-CG breakdown test (`so_bicg.cc`, `|r1.r2|`,
  `|d2^T A d1| < 1e-16`) and the convergence floor (`bicg_error_minimum`
  1e-16, `eps_stop` 1e-16 on `|r|^2`) are absolute, so both depend on the
  scale of the right-hand side. On tiny systems the false breakdowns are
  rescued by the exact direct band LU retry (C74), at the cost of a band
  factorization per solve.
- **Measured** (scratch variant, 2026-10-13): a scale-relative test
  (`|r1.r2| < 1e-14 |r1||r2|`, `|dAd| < 1e-14 |d2||Ad1|`) removes every
  breakdown on validation_1 (600/600 solves converge, result equal to
  1e-7, -10 % wall) with no corpus classification change (956 runs), but
  the converged Bi-CG then stops on the absolute floor and is LESS
  accurate than the LU retry: examp16 `node_dof` vs SuperLU 0.43 % ->
  2.45 % (`post_line_dof` sigyy 5.3e-4 -> 5.8e-3); validation_1 with the
  inflow scaled by 1e-5, velx at t = 0.05: HEAD 9.77556e-7 = SuperLU
  9.77517e-7, variant 8.96098e-7 (8 %). A restart-on-breakdown variant
  failed (14671 restarts, 200 solves still broke down). The Professional
  Bi-CG shows the same floor behaviour (8.9888e-7 vs PARDISO
  9.77518e-7); its default solver is PARDISO.
- **Next step**: scale-relative breakdown test TOGETHER with a relative
  convergence floor (e.g. `||r||/||b||`, or `bicg_error_minimum` scaled
  by `|b|^2`), gated by a corpus A/B against a direct solver (SuperLU)
  on the 38 inputs that break down today; check the cost of the extra
  Bi-CG iterations against the saved LU factorizations.
- **Owner notes** (2026-10-13): the owner suggested a stagnation stop
  (relative improvement between iterations). Planned design: main
  criterion `||r||/||b|| < tol`; stagnation detection over a window of k
  iterations on the BEST residual, returning the best iterate (the Bi-CG
  residual is non-monotone); a floor of about `eps*||A||*||x||` instead
  of the absolute 1e-16; and a true-residual `b - Ax` check at the stop
  (the recursive residual drifts near machine precision).
- **Resolution** (owner decision 2026-10-15, fix `d567672`): main
  criterion `|r|^2 < bicg_error |r0|^2` (relative); the absolute
  `bicg_error_minimum` only when the user gives the record; round-off
  floor `100 eps |A| |x|`; scale-relative Bi-CG breakdown (1e-14); a floor
  or breakdown stop is kept when `|b-Ax|/|b| <= 1e-6`, otherwise SuperLU,
  then the band LU within a memory budget (clear error instead of
  bad_alloc). No stagnation window and no best-iterate return (prototype:
  a 10-iteration window raised the retries 26 -> 10958). The CG takes the
  same stop (measured safe: E15 1.4e-6 -> 4.0e-11 off SuperLU, +7 % wall).
  Evidence: 86 Bi-CG corpus inputs, retries 7378 -> 13, 19 inputs 10x
  closer to SuperLU, 0 farther; large meshes C23 91.9 -> 18.6 s, C36
  bad_alloc -> 83 s; omen A/B 965 runs, 63 of the 84 moved runs toward
  SuperLU (details: dev/solve_iterative_bicg, SEGUIMIENTO). Test
  `bicg_relative_stop1`. DIVERGENCES C77 -> class A.
- **Pending owner choice** (item 2 of the decision): BiCGStab measured
  9-13 % less solver time on the large meshes but 71 SuperLU retries
  (29 max_iter solves) and 2 inputs farther from SuperLU on the 86-input
  corpus (`tquad6_iface` 5e-8 -> 7.5e-3); GMRES(30) clearly worse (379
  max_iter solves, 4.8x products). The classical Bi-CG with the relative
  stop is kept until the owner decides.
- **Severity**: LOW. **Effort**: medium.

### C81. ILU / AMG preconditioning of the iterative solver (F6, owner 2026-10-15) — OPEN (F6)

- **Context**: the CG / Bi-CG use only the symmetric Jacobi scaling
  `P = diag(A)^-1/2`. With the C77 relative stop the iteration counts
  grew (+66 % Bi-CG iterations on the 86 Bi-CG corpus inputs, CG E15 127
  -> 295).
- **Proposal**: an incomplete LU (ILU(0) / ILUT) or algebraic multigrid
  preconditioner on the assembled CSR operator (`so_bicg.cc`
  `bicg_csr_build()`), measured in the F6 profiling phase; it may also
  make BiCGStab robust enough to replace Bi-CG (C77 pending choice).
- **BiCGStab notes for the revisit** (C77 measurement, 2026-10-15): with
  Jacobi scaling only, BiCGStab stagnated on `distri3`, `bimet1` and
  `beam2d_2` and spent the full `max_iter = 10*n` before the direct
  retry; detect stagnation earlier so a hopeless solve reaches the
  direct retry sooner. On `tquad6_iface` it stopped at the round-off
  floor with a true residual below `BICG_FLOOR_ACCEPT` (1e-6) and the
  accepted iterate moved the result from 5e-8 to 7.5e-3 off SuperLU:
  review whether 1e-6 is too loose for interface systems.
- **Owner rationale** (2026-10-15): in large meshes with many steps even
  small per-solve gains compound.
- **Severity**: LOW (performance). **Effort**: medium/large.

### C73. Convection of the integration-point state of history materials (F4 T6f C72, 2026-10-13) — OPEN (optional improvement)

- **Context**: on a mesh fixed in space with `convection_apply -yes`
  neither the GNU nor the Professional 02-08-2026 convects the
  integration-point state (stress, plastic strains, kappa, damage,
  hypoplastic / Cam clay history, Maxwell stresses): it stays frozen at
  the integration points while the material flows through the mesh
  (probes, Engram 2775). For a fluid this is harmless; for a solid with
  memory (Eulerian forming, steady rolling) it is wrong physics.
- **Proposal**: convect the integration-point state (for example by
  re-interpolating it upstream along the mesh velocity each step, or by
  convecting a recovered nodal field and interpolating it back). This
  would be a deliberate class A divergence (GNU better than the
  Professional) and needs a project-owner decision.
- **Severity**: LOW (no corpus model relies on Eulerian solid history).
  **Effort**: large.

### C72. Convected (Eulerian) stress path (F4 T6f K21 WU1b, 2026-10-12; nodal output resolved F4 T6f C72, 2026-10-13) — OPEN (Professional-side residuals)

- **Corrected premise (2026-10-13, Engram 2775)**: the Professional
  02-08-2026 does NOT convect the integration-point stress. A sigzz,
  sigxy or isotropic-pressure sin^2 pulse (node_dof +
  `global_element_dof_from_node_dof -yes`) stays bit-frozen at the
  integration points and at the nodes for t = 0..2 (elastic and fluid,
  Stokes `-yes`/`-no`, velocity prescribed or free, uniform or parabolic
  shear, `-updated` or `-updated_linear`). The earlier wording ("the
  Professional convects / sweeps the stress") was wrong.
- **Resolved (`c1a1ae9`)**: the GNU nodal stress output. With the stress
  dof convected (`stress_recovery_absolute` 0) the lumped mass of the
  incremental update used the per-element recovery weight F/volume
  (`miscel.cc` `sri_stress_recovery_weight()`, `general.cc` lumped
  block) while the convective term used h, so the nodal stress moved at
  about v*volume_ip (0.055 for volume_ip 0.0625 and 0.011 for 0.0125 at
  v = 1). Fix: the convected stress dof is recovered absolutely, as
  everywhere else (the Professional's node stress is the recovery of the
  integration-point stress); the convective terms are skipped for it.
  `test-2026/convected_stress_recovery1` (40 quad4, uniform v = 1,
  elastic, sigzz pulse): post_point 0.5 / 1.0 / 0.5 at x = 1.5 / 2.0 /
  2.5, t = 2 = Professional; before 0.348732 / 0.904484 / 0.6481. The
  velocities never depend on the nodal stress (A/B on every input with
  `convection_apply -yes`: velocities bit-identical, nodal stresses move
  up to 1.3 % of their maximum on `validation_1` / `stokes_*`). The
  incremental path remains only with an independent strain unknown
  (`materi_strain_total`, e.g. sfnet `strconv1`, unchanged byte for byte).
- **Remaining (Professional side; (1) class D, (2)-(3) class B)**: (1) `validation_1` drift
  of the Professional (-0.0391 / -0.0179 / +0.0034 at t = 1 / 2 / 3,
  dt-converged), independent of `global_element_dof_apply`,
  stabilization, predictor and density; the GNU is steady at -0.057,
  inside the GAMM range -0.066 .. -0.018 (DIVERGENCES C72-V1, owner
  decision 2026-10-13: Pro-suspect; class D since 2026-10-18, decided,
  no longer open; guard `test-2026/validation1_steady1`). (2)
  Stokes flow: the Professional `validation_1` without the stokes record
  0.1465713 vs GNU 0.144159 = Pro `convection_apply -no` 0.1441588, a
  convection-gated, density-independent residual
  (`test-2026/stokes_default1` tolerates it, 0.003). (3) Original
  test-2014 `cd_method` rigid rotation (`-ra 1 4 -ra` node list): Pro
  sigxy 10.78, 9.89, 6.99, 4.22, 1.25, then 0 from t = 0.6; GNU 0.744 /
  3.642 at t = 1, not converged in dt or mesh. Both codes accumulate the
  compressibility pressure at the integration points linearly in time
  (not convected). Mechanisms of (1)-(3) not identified.
- **Severity**: LOW. **Effort**: medium (diagnosing a Professional-side
  effect without its source).

### C80.1. Selective / B-bar integration or an adaptive reduced factor (owner 2026-10-16, F4 T6f C80) — RESOLVED 2026-10-16 (F4 T6f C80.1)

- **Resolution** (F4 T6f C80.1, 2026-10-16, `b47c39c`): owner decision 2026-10-16:
  extend the GNU-only `group_element_selective_reduced_integration` with
  `-volumetric` (B-bar: volumetric part κ m m^T at the 1-point centroid,
  deviatoric part full Gauss, B-bar stored stress) and `-both` (shear +
  volumetric); `-yes` unchanged; opt-in only. Analytic: Lamé ν 0.4999
  8x16 −0.19 % (Gauss −80 %, Pro f=10 −11 %), no pressure checkerboard,
  no hourglass. Soft on 1-element-thick bending meshes (exact element
  behaviour, documented). Not done (owner questions): plasticity (scope
  stays linear elastic), axisymmetry (refused), the λ/μ split variant
  (better coarse bending, measured), an adaptive reduced factor.
  Tests `test-2026/sri_*` (4). Docs: user/dev
  `group_element_selective_reduced_integration`.
- **C80.2 follow-up** (F4 T6f C80.2, 2026-10-16, `5993499`, still closed): owner
  decision 2026-10-16: the λ/μ split added as `-lambda` (Hughes 1977:
  only λ·m·mᵀ at the centroid point, the 2μ part full Gauss; isotropic
  linear elastic only, non-isotropic groups refused, plane stress
  dropped as for `-volumetric`), `-both` removed (input error; enum
  slot reused by LAMBDA, numbering unchanged). Gate (owner): measured
  on the C80.1 problems before shipping, passed: Lamé ν 0.4999 8x16
  −0.24 % (−0.19 % `-volumetric`), no element-level checkerboard
  (0.1 %), cantilever 10x1 ν=0.3 −4.9 % (`-volumetric` +22 %), no
  hourglass, rigid-only zero eigenvalues (quad4 3, hex8 6). Tests
  `test-2026/sri_lambda_*` (3), `sri_both_removed_error1`.

- **Context**: the C80 `-reduced` rule (= Professional) applies one
  factor to the whole stress. The analytic study of C80 (plane-strain
  cantilever, thick cylinder, point load; tip / inner-radius error vs
  exact): f = 1 is stiffer than Gauss in bending (20x2 -46 % vs -11 %)
  and locks volumetrically (nu 0.49, 8x16: -11 %); f = 10 (the manual
  default) overshoots bending (20x2 +31 %, 40x4 +6.3 %, one element
  thick 10x1 +2100 % hourglass, zig-zag under a point load, no solver
  warning) but removes the volumetric locking (nu 0.49 -0.13 %, nu
  0.4999 -11 % vs Gauss -80 %); f = 2 is best in bending (40x4 -0.7 %,
  20x2 -2.8 %) and still locks near nu 0.5 (0.4999: -75 %).
- **Candidate**: selective / B-bar integration (volumetric part at 1
  point, deviatoric part at f ~ 2 or Gauss), or an adaptive f(nu). An
  explicit per-group opt-in, NEVER the default. It is a GNU extension
  (class C): a new keyword or value, to be decided by the owner before
  the work starts. Relation with the existing
  `group_element_selective_reduced_integration` (shear-only SRI) to be
  decided too.
- **Verification**: the C80 analytic cases (scripts in the C80 work
  notes) plus patch tests; no Professional reference.
- **Severity**: LOW. **Effort**: medium.

### C80. `group_integration_method -reduced` refused by the GNU (F4 T6f C79, 2026-10-13) — RESOLVED 2026-10-16 (F4 T6f C80)

- **Resolution** (F4 T6f C80, 2026-10-16, `995de49`): implemented after
  135 Professional 02-08-2026 probes (residual <= 2e-13). quad4/hex8 2x2
  (x2): points at iso +-1/f, equal weights; `-minimal`: 1 point. bar2
  default 1 point, `-maximal` Lobatto; bar3/quad9/hex27 Lobatto -1,0,1
  (factor ignored), `-minimal` 2-point Gauss; tria3/tet4/prism6 ignore
  the method. Nodal recovery of the quad/hex rules: physical linear fit
  plus the arithmetic mean, also with the points at the nodes (unlike
  `-lobatto`). Axisymmetric: an element with a node on the axis is
  integrated with Gauss for any factor (Pro probe at f = 1 and f = 2).
  Default factor: the GNU applies the manual 10, the Professional 1
  (owner 2026-10-16: deliberate divergence class A, DIVERGENCES C80).
  Tests `test-2026/reduced_quad4_f2_1`, `reduced_quad4_f10_1`,
  `reduced_hex8_f2_1`, `reduced_quad9_minimal1`, `reduced_quad9_maximal1`,
  `reduced_quad4_axisym1` (Pro targets 4e-6), `reduced_default_factor1`
  (GNU default): RED rc 1, GREEN, every node_dof = Pro within 5e-12. Not
  done: bar3 `-reduced` recovery (Pro linear fit, GNU 1.3-2.3 % off; C67, resolved 2026-10-16);
  the Professional element-average fallback at very large factor (1e6),
  not investigated. Follow-up: C80.1 (selective / B-bar, opt-in).

- **Symptom**: the Professional manual (6.623) lists `-reduced` (a
  weighted combination of the Lobatto and the 1-point rule, with a
  factor record) and the Professional 02-08-2026 runs a quad4 group with
  it; the GNU stops with "I do not know : -reduced". The GNU user manual
  page `integration_quadrature_and_stress_recovery` listed `-reduced` as
  accepted (corrected 2026-10-13).
- **Fix candidate**: implement the reduced rule after measuring the
  Professional stiffness and recovery with it, or document it as not
  supported.
- **Severity**: LOW (no corpus input uses it). **Effort**: medium.

### C79. Nodal recovery of the Lobatto quad4 is volume weighted; the Professional takes the arithmetic mean (F4 T6f C71, 2026-10-13) — RESOLVED 2026-10-13 (F4 T6f C79)

- **Resolution** (F4 T6f C79, 2026-10-13): Professional 02-08-2026
  probes (prescribed nonlinear velocity at every node, distorted meshes
  of unequal element size, `element_intpnt_dof` -> `node_dof`, residual
  <= 5e-12 of a field scale 13-19): for the Lobatto rules with the points
  at the nodes (quad4 2x2, quad9 3x3, hex8 2x2x2, hex27 3x3x3; plane and
  axisymmetric off the axis) the per-element value is the point value at
  the node (the physical linear fit is 0.39-0.46 off) and the node is
  the arithmetic mean; for the quad9 2x2 / hex27 2x2x2 Lobatto
  `-minimal` (corner points) it is the physical linear fit `1,x,y(,z)`,
  as for the 2x2 Gauss rule (the point value is 0.18-0.19 off, also on
  a rectangle); quad4 / hex8 1 point already matched. GNU cause:
  `sri_stress_recovery_weight()` returned `h_inol` for every Lobatto rule
  (times the point volume in the lumped equation = volume weighting), and
  for the quad9/hex27 `-minimal` rules `h = 0` at the non-corner nodes, so
  the `materi.cc` stress block skipped them (sigxx 0 instead of up to
  16.5). Never measured before (`06e175d` kept `h` for the Lobatto rules);
  no earlier test or corpus target used `-lobatto` (sfnet `tube1` was
  NaN until C71). Fix (`miscel.cc` `recovery_fit_kind()` kind 3 point
  value, Lobatto point positions in `recovery_fit_shape()`; `polynom.cc`
  `pol()` sets `stress_recovery_gauss` for the `-minimal` Lobatto rules
  and a new `stress_recovery_lobatto`, both thread local). Tests
  `test-2026/recovery_lobatto_quad4_mean1`,
  `recovery_lobatto_quad4_axisym_mean1`, `recovery_lobatto_hex8_mean1`,
  `recovery_lobatto_quad9_minimal1`, `recovery_lobatto_hex27_minimal1`
  (Professional targets 4e-6; RED rc 1, GREEN). sfnet `tube1` (modern
  keywords): every node_dof = Professional within 6e-8 (before: the
  r=1.25 column, node (1.25,0) sigxx -494.42 vs -407.64). Interface
  elements are not involved (`interface.cc` assembles forces only). The
  C71 leftover (`calcul_force.cc` `msf_element_rule()` replica of the
  on-axis switch) was closed in the same unit: test
  `test-2026/msf_axisym_lobatto_axis1` (momx_sig 1/6 = the `-gauss` run;
  before 0.5). Still unmeasured Lobatto rules: bar2/bar3, quad16, hex64
  (C67; bars and quad16 resolved 2026-10-16, hex64 rejected by both codes). Found: `group_integration_method -reduced` refused by the GNU
  (C80).

- **Symptom**: with `group_integration_method -lobatto` the GNU nodal
  stress of a quad4 is the volume-weighted average of the adjacent
  elements' point values at the node (`sri_stress_recovery_weight()`
  returns `h` for the node-containing rules, and the lumped sigma-dof
  equation multiplies it by the point volume). The Professional
  02-08-2026 takes the arithmetic mean of the elements' values at the
  node, as for its Gauss rules. Two plane quad4 of widths 1 and 2 sharing
  an edge (quadratic prescribed velocity): shared node sigxx GNU 6.525
  (weights 1/3, 2/3), Professional 6.15 = mean(5.025, 7.275); the same
  for axisymmetric off-axis elements (7.9 vs 7.525). Equal sizes agree.
  quad9 Lobatto: equal in the probe (the element values at the shared
  nodes coincide), so not discriminated.
- **Visible in**: sfnet `tube1` (`-lobatto` + axisymmetric) after C71:
  velocities = Professional (7e-9), nodal stresses = Professional except
  the r=1.25 column, where a Gauss axis element (weight 1, C71) meets a
  Lobatto element (weight = point volume 2*PI*r*w*detJ ~ 15): node
  (1.25,0) sigxx -494.42 vs -407.64 = mean(-308.77, -506.52).
- **Fix candidate**: return `h/volume` (each element weight 1) for the
  quad4 Lobatto rule, after measuring hex8 / quad9 / bar Lobatto and the
  interface elements (which also integrate at the nodes).
- **Severity**: LOW (only explicit `-lobatto` groups with unequal
  element volumes at a node). **Effort**: cheap fix, medium measurement.

### C71. sfnet `tube1` (`-lobatto` + axisymmetric) gives NaN in the GNU (F4 T6f K21, 2026-10-12) — RESOLVED 2026-10-13 (F4 T6f C71)

- **Resolution** (F4 T6f C71, 2026-10-13): first NaN at `polynom.cc`
  `pol()`: the Lobatto points of an element touching the axis lie at
  r=0, so the hoop entry `new_b += h/radius` is infinite (or 0/0) and the
  point volume `2*PI*r*w*detj` is 0; `B^T sig dV` = inf*0 = NaN in the
  velx right-hand side (`node_rhside` velx -nan on a one-element probe),
  then NaN velocities. Professional 02-08-2026 probes (one on-axis and
  one off-axis quad4/quad9 in a `-lobatto` axisymmetric group,
  `element_intpnt_dof`): the element touching the axis gets the Gauss
  points and values (identical to a `-gauss` run), the off-axis element
  keeps the Lobatto points. Switch on the node radius, absolute: 9e-10
  Gauss, 1e-9 Lobatto (scale 1e-6 with r=1e-10: Gauss; scale 1000 with
  r=1e-7: Lobatto), i.e. `|r| < EPS_COORD` (1e-9). The Professional
  manual (6.621, 6.623) does not mention it. Fix: `pol()` rebuilds the
  rule with Gauss when an axisymmetric Lobatto point has `|r| <
  EPS_COORD` (also sets `stress_recovery_gauss`, so the on-axis recovery
  rule of K21 applies). Test `test-2026/axisym_lobatto_axis1` (RED rc 1
  on 58af21e: node 3 vely 1.1, node 1 stresses 0; GREEN = Professional
  1e-9). sfnet `tube1` (modern keywords, `-updated_linear`): velocities
  = Professional within 7e-9; the remaining nodal stress difference at
  the r=1.25 column is the Lobatto recovery weighting, C79. The
  `post_calcul -materi_stress -force` replica of the rule
  (`calcul_force.cc` `msf_element_rule()`) applies the switch since
  F4 T6f C79 (2026-10-13).

- **Corpus**: sfnet `tube1` (`group_integration_method -lobatto`,
  `group_axisymmetric -yes`, h-refined quad4, viscous fluid).
- **Symptom**: GNU `node_dof` velocities NaN ("the iterative solver did
  not converge"); the Professional 02-08-2026 runs it (node 3 vely
  999.98, sigxx -158.25). The Pro `tube1` (default Gauss rule) runs in
  both programs and is a K21 test (`recovery_tube1_axis1`).
- **Cause candidate**: the Lobatto points on the axis have radius 0, so
  their axisymmetric volume `2*PI*r*w*detJ` is 0 (and terms divided by r
  are singular); how the Professional treats on-axis Lobatto points is
  not measured.
- **Next step**: one-element on-axis Lobatto probe in both programs
  (element_intpnt_dof, volumes), then decide.
- **Severity**: LOW. **Effort**: cheap probe, medium fix.

### C70. 3x3 quad9 and 3D recovery tables fit in natural, the Professional in physical coordinates (F4 T6f K21, 2026-10-12) — RESOLVED 2026-10-13 (F4 T6f C70)

- **Resolution** (F4 T6f C70, 2026-10-13): Professional 02-08-2026 probes
  (distorted / curved elements, `element_intpnt_dof` +
  `element_intpnt_coord` -> `element_dof`): quad9 3x3, hex8 2x2x2,
  hex27 3x3x3 and 2x2x2, tet10 5-point and prism6 6-point = linear least
  squares `1,x,y(,z)` in PHYSICAL coordinates (residual <= 1e-12 of the
  scale); hex8 1 point = the point value; axisymmetric quad9 (3x3 and
  2x2) with a node at r=0 = `sum h sig / sum h` (the GNU 2x2 used the
  linear fit there); node = arithmetic mean (1e-13). `miscel.cc`
  `recovery_fit_kind()` / `recovery_fit_shape()` / `recovery_fit_build()`
  generalised to 3D; the static natural table (F4 T9) removed. Tests
  `test-2026/recovery_hex8_distorted1`, `recovery_hex8_point_mean1`,
  `recovery_hex27_distorted1`, `recovery_tet10_curved1`,
  `recovery_prism6_linear1`, `recovery_quad9_curved1`,
  `recovery_quad9_axis1` (RED rc 1 on 421a3f4, GREEN). Side effect:
  `force11` hex27 nodal stresses = Professional (C57 recovery part).
  Found: the prism6 integration points differ (C75). Unmeasured rules
  stay in C67 (prism15, Lobatto, ...).


- **Symptom**: K21 WU1 probes (prescribed quadratic velocity at every
  node, viscous fluid, Professional `element_dof` vs the operators
  applied to its `element_intpnt_dof`): the Professional per-element
  value of a CURVED/distorted quad9 with 3x3 points and of a DISTORTED
  hex8 with 2x2x2 points is the uniform linear least squares in
  PHYSICAL coordinates (residual 1e-12), not the natural-coordinate fit
  of the GNU table (`sri_recovery_factor()` types 1 and 2): the natural
  fit is off by ~1.0 of a field scale 8 (quad9) and 1.1-1.5 of a scale 10
  (hex8). Equal on parallelograms / parallelepipeds (affine invariance),
  which is why the earlier probes matched.
- **Not changed in K21 WU1** (scope: 2D 2x2 / 1-point quad rules). The
  hex27 and tet10 tables are probably the same case (unmeasured; C57 hex27
  `force11` up to 3.5 % may be this). The geometry-dependent builder
  `recovery_fit_build()` (`miscel.cc`) generalises to 9 points directly;
  3D needs the 4-term fit and the 3D shape functions.
- **Next step**: extend `recovery_fit_kind()` to quad9 9-point (one line +
  test), then a 3D builder for hex8/hex27 (+ tet10 probe), corpus A/B.
- **Severity**: LOW (nodal output only). **Effort**: cheap (quad9) /
  medium (3D). (Original entry, kept for history.)

### C69. Plane-stress membrane: the yy/zz split of the plastic strain is first order in dt (F4 T6f, 2026-10-12) — DEFERRED (class B, F4 T6f LOW queue, owner plan 2026-10-17)

- **Deferred (2026-10-17)**: a plane-stress-consistent return plus a Professional eppyy/eppzz
  comparison (medium); GNU-only residual that converges with dt.

- **Context**: found while writing the guard tests of the deliberate
  divergences MEMBRANE-EPTZZ and MC-VERTEX-ZZ (owner decision 2026-10-12,
  class A in `DIVERGENCES-PROFESSIONAL.md`). Those two K15 items
  (`validation_3` eptzz, `mmc_tension` sigzz/eppzz) are closed as
  deliberate divergences with test-2026 `membrane_eptzz_total1/2` and
  `mmc_tension_vertex1/2`; this item is the GNU-only residual.
- **Symptom / evidence**: single quad4, von Mises 24.3 + kap hardening
  (slope 224), E 7000, nu 0.2, `group_materi_membrane -yes`, uniaxial
  stretch eptxx 0.02. Analytic (uniaxial Prandtl-Reuss): eppyy = eppzz =
  -0.0080414, eptzz -0.0088249. GNU eppyy / eppzz / eptzz:
  dt 0.1 -0.009654 / -0.006426 / -0.0072102;
  dt 0.01 -0.008236 / -0.007847 / -0.0086300;
  dt 0.001 -0.008061 / -0.008022 / -0.0088050. eppxx 0.01608, sigxx
  27.42, sigzz ~0 and tr(epp) = 0 at every dt; eptzz = -nu sigxx/E +
  eppzz holds (the stored total is consistent).
- **Root cause (probable)**: `membrane_apply()` (`membrane.cc`, called in
  `stress.cc` `set_stress()` ~1871) makes sigzz zero with an ELASTIC zz
  strain correction on the trial stress, before the plastic return; the
  return of a step therefore does not see an exact plane-stress state.
- **Fix candidate**: a plane-stress-consistent return (iterate the zz
  correction with the plastic return, or the reduced plane-stress return
  mapping). Needs a Professional comparison of eppyy/eppzz first.
- **Severity**: LOW (converges with dt; the total eptzz is stored).
  **Effort**: medium.

### C68. Other integration-point state variables read back from the nodes (F4 T6f WU2, 2026-10-11) — RESOLVED 2026-10-12

- **Symptom / evidence**: the ELEMENT_DOF masks of `elem.cc` `elem()`
  (read ~636-650, write ~812-825) restore/store only the history
  variables, stresses, the epe/epp/ept/epi strains, svloc/svnonloc and
  (since F4 T6f WU2) kappa. Every other state that the constitutive
  update reads from `old_unknowns` is therefore the h-lumped NODAL value
  interpolated back to the point, the defect fixed for kappa in WU2:
  `materi.cc` `old_kapsh` (`kapsh_indx`, ~287), `old_cap1pc`
  (`cap1_indx`, ~292), `old_rho` (`rho_indx`, ~265, kinematic
  hardening), `old_damage` (`dam_indx`, ~297) and the Maxwell stresses
  `old_msig` (`mstres_indx`, ~262). `old_f` and `old_substeps` are only
  used for their own nodal output (not state inputs).
- **Not changed**: no Professional evidence per variable yet (kappa was
  validated on `validation_3`). Candidate tests: Pro `other` models with
  `materi_plasti_kappa_shear` (MC hardening-softening), cap1, rho,
  damage, Maxwell.
- **Fix candidate**: the same two-line mask extension + the IP write in
  `materi.cc`, one variable at a time with a two-element probe vs the
  Professional and a corpus A/B.
- **Severity**: MEDIUM (wrong hardening history on non-uniform fields).
  **Effort**: cheap per variable, medium with the evidence.
- **Status**: **RESOLVED 2026-10-12 (F4 T6f, `70235c6`)**: the manual
  (4.10, 4.15, 4.17, 4.25, 4.27) treats the five as point history, like
  kappa (4.24), so it is a GNU defect fix. `elem.cc` masks (read and
  write) gain kapsh, cap1, dam, the rho range and the Maxwell range
  (feature-guarded); `materi.cc` the old h-weighted blocks are replaced
  by one block after the kappa block, outside the h guard, with
  `strain_recovery_weight()` and the IP write; `general.cc` denominator
  weight. Verified on four-tria3 probes against the Professional
  02-08-2026 (tests `test-2026/damage_ip_state1`, `rho_ip_state1`,
  `kappa_shear_ip_state1`, `maxwell_ip_state1`, all RED rc 1 on the
  previous binary, GREEN rc 0): damage reactions 1.9e-1 -> 5e-13
  (relative), rho 1.1e-1 -> 7e-6, kapsh 1.3e-2 -> 1.5e-4, Maxwell nodal
  msigyy 8.4e-2 -> 5e-4. cap1 is changed for consistency with the
  manual only (no Professional evidence, C68b). Corpus A/B on omen
  (`f4t6f_c68_runs`): no rc or class change, no node_dof change above
  1e-10. Residuals found: C68a, C68b, C68c, C68d.

### C68a. Professional Maxwell chain differs from the analytic solution (F4 T6f C68, 2026-10-12) — DECIDED 2026-10-18 (class D, Professional-suspect)

- **Decision (owner 2026-10-18)**: class D, GNU kept (DIVERGENCES C68a,
  SEGUIMIENTO Descarte). Analytic: dmsig/dt = E_c r - msig/tau, msig(0) =
  0, so msig(t) = E_c tau r (1 - exp(-t/tau)) = 7e3 x 0.3 x (-0.01) x
  (1 - exp(-1/0.3)) = -20.2508461397; GNU equals it to 12 digits at dt
  0.1 and 0.01. Professional 02-08-2026: -23.8131 (dt 0.1), -20.5902
  (dt 0.01), -20.2846 (dt 0.001), a first-order time-discretisation
  error; at dt 0.1 it exceeds the bound E_c tau |r| = 21. Guard test
  `test-2026/maxwell_relaxation_analytic1` (1e-6 relative; a wrong
  target gives rc 1). Probes: scratchpad `c68a/`.

- **Symptom / evidence**: homogeneous one-chain relaxation probe (tau
  0.3, constant strain rate): Professional 02-08-2026 msig -23.81; the
  analytic E*tau*rate*(1-exp(-t/tau)) and GNU -20.25. On a
  non-homogeneous field the Professional also gives msigxx != 0 where
  GNU gives ~0. The C68 test `maxwell_ip_state1` matches the
  Professional within 5e-4 nonetheless (the difference is the chain
  law, not the IP state).
- **Classification**: class D since 2026-10-18 (was class B,
  Professional-suspect, not decided).
- **Severity**: LOW. **Effort**: done (guard test + class D row).

### C68b. cap1 model: GNU and Professional differ at model level (F4 T6f C68, 2026-10-12) — DEFERRED (class B, F4 T6f LOW queue, owner plan 2026-10-17)

- **Deferred (2026-10-17)**: model-level cap1 investigation (medium/expensive), no corpus
  case affected.

- **Symptom / evidence**: cap1 probe: Professional c1his -17 vs GNU pc
  108; reactions differ by a factor 6. The cap1 IP-state change of C68
  could therefore not be verified against the Professional.
- **Classification**: class B, cause unknown.
- **Severity**: LOW. **Effort**: medium/expensive (cap1 law and its
  history variable, `group_materi_plasti_cap1`).

### C68c. Mohr-Coulomb hardening-softening on kappa_shear: 1.2x off the Professional (F4 T6f C68, 2026-10-12) — DEFERRED (class B, F4 T6f LOW queue, owner plan 2026-10-17)

- **Deferred (2026-10-17)**: cause unknown, homogeneous model-level difference (medium).

- **Symptom / evidence**: homogeneous control `c68_kapsh2*` (scratch
  probes of the C68 investigation): GNU and the Professional differ by a
  factor 1.2 with `group_materi_plasti_mohr_coul_hardening_softening`
  driven by `-kapsh`; homogeneous, so not the IP-state defect (the von
  Mises kapsh test `kappa_shear_ip_state1` matches within 1.5e-4).
- **Classification**: class B, cause unknown.
- **Severity**: LOW. **Effort**: medium.

### C68d. `materi.cc` old-state read uses the Maxwell loop index (F4 T6f C68, 2026-10-12) — RESOLVED 2026-10-12 (`2327083`)

- **Symptom**: in the old-state read loop of `materi()` the Maxwell loop
  overwrites `indx` with `m*MDIM*MDIM + idim*MDIM + jdim`, and the
  following `old_rho[indx]`, `old_epi[indx]` and `old_epp[indx]` reads
  use it. With one chain the value equals the intended index; with two
  or more chains and rho, epi or epp initialised they write past the
  9-entry arrays. Found by code reading; untested combination, no
  Professional comparison involved.
- **Fix candidate**: recompute `indx = idim*MDIM + jdim` after the
  Maxwell loop (or use a separate index) plus a two-chain + plasticity
  test.
- **Severity**: LOW (latent; memory safety when triggered).
  **Effort**: cheap.
- **Resolution** (F4 T6f C68d, `2327083`): confirmed by code reading
  (only the old-state read loop had the pattern; the C68 write block
  and the rotation loops use explicit indices) and by a test. With two
  chains the rho/epi/epp reads went past their 9-entry slices into
  `new_rho`, `new_epi` and `inc_epp` (inside the single `dbl_array`
  allocation, so an ASan/UBSan build reports nothing) and the real
  slots stayed 0: rho and epp never accumulated. The Maxwell read now
  indexes `old_msig[m*MDIM*MDIM + indx]`. Test `test-2026/
  maxwell_two_chains_state1` (von Mises + kinematic hardening, two
  identical chains = one chain with the summed stiffness, exact in the
  Professional 02-08-2026 too): previous binary rc 1 (node 1 reaction
  -13.32 vs -19.17), fix rc 0 (1e-11). The remaining GNU-vs-Professional
  difference of that model (one chain: reaction 2.5e-3, rhoyy 7e-3,
  msigyy 9 % relative) is the Maxwell relaxation item C68a.

### C67. Nodal recovery of the unmeasured element rules (F4 T6f, 2026-10-11) — RESOLVED 2026-10-16 (F4 T6f C67 `d875874` and C67 remainder `dc4a59a`)

- **Symptom**: F4 T6f WU1 aligned the nodal recovery of tria3 (1 point,
  3-point nodal rule), tria6, tet4 (1 point) and bar2 (1 point) with the
  Professional arithmetic mean (two-element probes, exact). The other
  rules outside the Gauss table still return `h_inol` (volume-weighted
  across elements): tet4 `-maximal`, bar2 with 2 points, bar3/bar4, the
  quad/hex Lobatto rules (`group_integration_method -lobatto`), tria3/tet4
  axisymmetric variants were not probed separately. Since C70
  (2026-10-13) hex8 1 point, hex27 2x2x2 and prism6 are measured and
  aligned; still `h_inol`: prism15, quad16/hex64, the Lobatto rules,
  tet4 `-maximal`, bar2 2 points, bar3/bar4. Since C79 (2026-10-13) the
  quad4/quad9/hex8/hex27 Lobatto rules are measured and aligned; the bar
  and quad16/hex64 Lobatto rules remain. C80 (2026-10-16): the
  Professional recovers the bar3 `-reduced` rules (Lobatto 3 points or
  Gauss 2 points) with the linear fit; the GNU keeps the `-lobatto` /
  Gauss bar recovery (1.3-2.3 % off on the C80 probe).
- **Bars and quad16 resolved** (F4 T6f C67, 2026-10-16, `d875874`): 65
  Professional probes (bar2/bar3/bar4 in 1D, quad16; Gauss, `-lobatto`,
  `-reduced` f=2/4, default/`-minimal`/`-maximal`; distorted, unequal
  sizes, every node). Rule = the quad/hex one: Gauss, `-reduced` and
  Lobatto `-minimal` = per-element physical linear fit `1,x(,y)`
  (residual 5e-13; natural fit 1e-2); Lobatto `-maximal` = point value
  (the bar4/quad16 interior nodes are at the Lobatto positions in both
  codes); arithmetic mean (length-weighted 1e-2 off); quad16 on the axis
  = `sum h sig / sum h` after the C71 Gauss switch. GNU cause: `h`
  weights (volume-weighted mean; bar3 Lobatto `-minimal` mid nodes 0),
  not intended (unmeasured fallback); no existing test covered these
  rules. Fix in `miscel.cc` (dimension-generic fit, nodes at the Lobatto
  positions) and `polynom.cc` (every `-reduced` Lobatto rule fitted);
  every probe node = Professional within 4e-12; 10 tests
  `test-2026/recovery_bar*`, `recovery_quad16_*`. `hex64` is rejected by
  both codes; truss/beam elements have no recovered nodal stress.
- **Remainder resolved** (2026-10-16, `dc4a59a`): enumeration of every
  rule `pol()` produces against the recovery dispatch; on the `h`
  fallback were tet4 `-maximal`, prism15 (every value) and prism6
  `-minimal` (C76). 81 Professional probes (tria3, tria6, tria3
  axisymmetric off/on the axis, tria6 axisymmetric, tet4, tet10, prism6,
  prism15; default/`-gauss`/`-lobatto` x default/`-minimal`/`-maximal`):
  tet4 `-maximal` = the 5-point Keast rule (manual 6.623; GNU had 4 nodal
  points) and prism15 = physical linear fit (5e-13), arithmetic mean.
  The prism15 7-point triangle table of `pol()` had (b,b)/(d,d) instead
  of (a,a)/(c,c) (a point outside the triangle, monomial error 0.025):
  corrected, elastic displacements 3.7 % -> 3e-11 from the
  Professional. tria3/tet4 axisymmetric: `pol()` forces 1 point
  (= Professional for the default; the `-maximal` difference is C83);
  tet4 never axisymmetric (3D). Every probe node of a changed rule =
  Professional within 6e-12. Tests `test-2026/recovery_tet4_maximal1`,
  `recovery_prism15_fit1`, `prism15_elastic_shear1`. Found: C66
  re-scoped (tria6 rule), C82, C83.
- **Fix candidate**: one two-element probe per rule (different size and
  stiffness), compare the Professional `node_dof` with the mean and the
  volume-weighted mean of its `element_dof`; extend
  `sri_stress_recovery_weight()` only for the rules proven to be the mean.
- **Severity**: LOW. **Effort**: cheap.

### C88. Spring element: elongation and `-spring1` right-hand side (F4 T6f C86, 2026-10-19) — RESOLVED 2026-10-19 (owner: fix)

- **Fix (F4 T6f C88)**: `spring()` (`spring.cc`) takes the elongation as
  the relative displacement projected on the direction (the record, or
  the INITIAL first -> second node axis), with the signed tangent
  `k dt dir dir^T` and the `-spring1` right-hand side / matrix at the
  right indices with a positive stiffness. Professional 02-08-2026
  probes (2D, k 10): oblique spring2 2.828427124746 (1 and 10 steps),
  rotated spring2 0 for every memory value, spring1 -2 with node_rhside
  +2, two oblique springs on a free node = analytic; the C86 probes now
  match (spring1 interior 100 / 5.496711 / -20, 3D oblique spring
  22.384 and section 100 / 10 / -4 / -40 / -6.3836). Tests
  `spring2_oblique_projection1`, `spring2_rotation_initial_axis1`,
  `spring1_projection1`, `spring_oblique_free1`,
  `spring2_updated_rotation1`, `post_element_force_2d_spring1_interior1`,
  `post_element_force_3d_truss_spring1`.

- **Measured** (Professional 02-08-2026 vs GNU, C86 probes,
  `post_element_force` strip): spring2 along its direction and along
  the motion: equal (55.385). Oblique spring2 (3D, direction 1 0 0, node
  6 -> 20): Professional 22.384, GNU -44.437 after one step, 8.093 after
  100 steps; aligned 3D spring2 6 -> 10 with bending: Professional
  -3.975, GNU -0.172; spring1 (direction 0 1, interior node): Professional
  4.503, GNU 21.168 (one step) / 6.256 (100 steps); spring1 (direction
  1 0, x loads only): Professional 1.887, GNU 9.899, and the GNU section
  force 40.2 instead of ~98 (no equilibrium).
- **Cause (code reading, not yet fixed)**: `spring.cc` takes the
  elongation as the change of the node distance (`array_distance`, the
  full motion, unsigned for spring1) while the force and the tangent act
  along `group_spring_direction` (`k dir dir^T`): the residual is not
  the derivative of the tangent, one step leaves an unbalanced force
  and the result is step-dependent. The Professional behaves as the
  projection `dir . (u2 - u1)` (linear, step-independent). In the
  `-spring1` branch `element_rhside[indx]` is written before `indx` is
  set for the current dof (stale index).
- **Severity**: MEDIUM (wrong spring forces for any spring that moves
  off its direction). **Effort**: cheap/medium (projection, the
  `-updated` memory keeps the rotating axis; spring1 index), plus corpus
  A/B.

### C89. Mohr-Coulomb interface dilatancy with `group_interface_gap` (F4 T6f C87, 2026-10-21) — DEFERRED (class B, owner 2026-11-01)

- **Measured** (Professional 02-08-2026, `iface_mc_dil` geometry: two
  quad4 blocks and a thick quad4 interface, c = phi = 0, psi 30 deg,
  right block pushed 0.1 and slid 5 over time 1): interface normal
  stress -1493.376 with gap 1.e20 or 10 for any step size; -2971.967
  with gap 0, 0.01 or 0.5 for any step size; 0 with gap -0.01 (a
  compressed interface reported open). The GNU after C87: equal for gap
  1.e20 / 10 and for one step with gap 0, 0.01 or 0.5; with two or more
  steps it gives -1493.376 for every gap below 10, and -2971.967 for one
  step with gap -0.01.
- **Cause (hypothesis)**: the Professional compares the gap with a
  normal measure that includes the dilatant opening (positive), not the
  accumulated effective strain. Not identified.
- **Tests**: `iface_mc_dil_1step` (default gap, one step, -2.97197 =
  Professional) and `iface_mc_dil` (explicit gap 1.e20, -1.49338 =
  Professional) are both pinned to the Professional.
- **Severity**: LOW. **Effort**: medium (probes of the per-step state).
- **Deferred** (owner 2026-11-01, closing F4 T6f): narrow scope (dilatancy
  plus a finite gap below 10 plus several steps) and a Professional
  answer that looks defective (0 stress on a compressed interface with
  gap -0.01); no manual rule decides it. Revisit after F5.

### C99. vtk and tecplot file names differ from the Professional (F5 W2 V2, 2026-11-02) — RESOLVED 2026-11-02

- **Measured** (Professional 02-08-2026, `control_post1` quad4 geometry):
  `control_print_vtk -yes` writes `<name>.vtk` (`-separate_index`:
  `<name>_000.vtk`), `control_print_tecplot` and `print_tecplot_calculation`
  append to `<name>.plt`. The GNU writes `tn<index>.vtk` and `tecplot.plt`
  (`print_vtk_calculation` already writes `<base>.vtk`). DIVERGENCES `C99`
  (class B).
- **Severity**: LOW. **Effort**: small (file names in `print_vt.cc` /
  `print_te.cc`; scripts and tests that read `tn*.vtk` / `tecplot.plt` must
  follow).
- **Resolution** (owner 2026-11-02 "sigamos a Pro", separate formats "el de
  Pro"): measured Professional names implemented in `print_vt.cc` /
  `print_te.cc`, including `-separate_index` / `-separate_sequential` (the
  GNU ignored them) and the truncation of `<name>.plt` at the start of a
  run. Against the manual 6.339 names: DIVERGENCES `VTK-TECPLOT-NAMES`
  (class A). Six test-2026 tests; the test-2014 `vtk_*1` inputs use
  `-separate_index`.

### C100. `global_element_dof_from_node_dof` default (F5 T2 V5b, 2026-11-03) — RESOLVED 2026-11-03

- **Measured** (Professional 02-08-2026, probes f5/v5b `d_*`): the manual
  (6.554) and the Professional default is `-no`: `node_dof` stresses of the
  input are not taken into the integration points. The GNU default stays
  `-yes` (its behaviour before the record existed). DIVERGENCES
  ELEMENT-DOF-FROM-NODE-DOF-DEFAULT (class B).
- **Severity**: LOW. **Effort**: small (initial value in `top.cc`), but
  about 200 suite inputs hold `node_dof` records and may change.
- **Status**: **RESOLVED 2026-11-03** (owner decision: "como manual y
  Pro"). Default `-no`. The first flip broke `hypo_masin_history_slots1`
  (initial state from `control_reset_dof`, no `node_dof` input), which the
  Professional runs: `data.cc` now writes the reset dofs into the
  integration points under `-no`, with the reset value evaluated at the
  point coordinates (= Professional element_intpnt_dof of
  control_reset_value_exponent1 / _logarithmic_first1). 24 suite inputs
  that give their initial state in `node_dof` (dof1, hypo1-4, cda1,
  smooth1, mhardsoil_*, m*diprisco*, mc_pressure_min*, mcap1, mcap1_comb,
  mstrain_cap, groundflow_total_pressure_tension,
  hardsoil_kappa_iterations1/2, hypo_wolff_igs_oedometer1/2) and
  convected_stress_recovery1 write
  `global_element_dof_from_node_dof -yes`; they fail the Professional
  parser for other keywords, so their Professional default run could not
  be compared. Test `global_element_dof_from_node_dof3` pins the default.

### C101. `global_node_dof_empty -yes` and the displacement (F5 T2 V5b, 2026-11-03) — RESOLVED 2026-11-03

- **Measured** (probes `c_def`, `c2_def`): the Professional zeroes every dof
  of a node only connected to empty elements except the displacement; the
  manual (6.552) says every dof. The GNU follows the Professional for now.
  DIVERGENCES GLOBAL-NODE-DOF-EMPTY-DISPLACEMENT (class B).
- **Severity**: LOW. **Effort**: small (`bounda.cc` `node_dof_empty_apply()`).
- **Status**: **RESOLVED 2026-11-03** (owner decision: keep the
  displacement, as the Professional and the current GNU; zeroing it risks,
  with an updated memory, nodes of a re-activated element lying inside
  other elements). No code change; DIVERGENCES row moved to class A.

### C102. Velocity of a node without stiffness (F5 T2 V5b, 2026-11-03)

- **Measured** (probes `c_no`, `c2_no`): the Professional gives a node that
  is only connected to empty elements zero velocity even with
  `global_node_dof_empty -no`; the GNU leaves it unchanged (node_lhside 0),
  so the displacement grows by v*dt. Pre-existing. DIVERGENCES
  NODE-NO-STIFFNESS-VELOCITY (class B).
- **Severity**: LOW. **Effort**: small/medium (solver rows with a zero diagonal).

### C103. `dof_element_dof` set (F5 T2 V5b, 2026-11-03)

- **Measured** (probes `g_w`, `g_f`, `g_vf`): the Professional marks every
  non-principal dof except the displacement `-yes` (work, f, void too); the
  GNU remembers only the constitutive states per element and marks those.
  DIVERGENCES DOF-ELEMENT-DOF-SET (class B).
- **Severity**: LOW. **Effort**: medium (per-element storage of the extra
  dofs in `elem.cc`, with a check of every reader).

### C104. `dof_limit` and prescribed values outside the limit (F5 T2 V5b, 2026-11-03)

- **Measured** (probe `b_m22`): the Professional element pass uses the
  unclamped prescribed values; the GNU post-solve material pass sees the
  clamped ones. Final `node_dof` equal. DIVERGENCES DOF-LIMIT-PRESCRIBED
  (class B).
- **Severity**: LOW (contradictory input). **Effort**: small.

### C105. Temperature slot of `element_intpnt_dof` (F5 T2 V5b, 2026-11-03)

- **Measured** (probe `b_m22`, also on c01bb20): with `condif_temperature`
  and mechanics the GNU writes 0 in the temperature slot of
  `element_intpnt_dof`; the Professional writes the interpolated
  temperature. Output only. DIVERGENCES ELEMENT-INTPNT-DOF-TEMPERATURE
  (class B).
- **Severity**: LOW. **Effort**: small/medium (find where the condif
  element path skips the integration-point PUT).

### C106. `control_reset_value_dof` combined with a spatial value record (C100 review, 2026-11-03)

- **Found** (code reading, `data.cc` `reset_value_at()`): the dof diagram
  and the spatial coefficients share the buffer `reset_value_diagram`; the
  spatial `db()` read of the first point overwrites the diagram, so the
  following points read the spatial coefficients as the dof diagram.
  Pre-existing (the node loop before C100 had the same order). No suite
  input combines both records; not compared with the Professional.
- **Severity**: LOW. **Effort**: small (separate buffer per record).

### C116. Two plasticity models in one group (tension + vonmises) differ from the Professional (F5 T3 D2, 2026-11-06)

- **Found:** 1D bar of `materi_plasti_visco_apply1` (E 2, strain 1.5) with `group_materi_plasti_tension 0 1.` and `group_materi_plasti_vonmises 0 0.8`: GNU sigxx 0.8283 / eppxx 1.0858, Professional 02-08-2026 sigxx 1 / eppxx 1; with `group_materi_plasti_visco_power 0 0.1 1.` GNU 2.6482 (= vonmises alone with visco), Professional 1.1987. Each model alone agrees (tension 1 / 1; vonmises 1.5333 / 1.5333; visco 2.6667 and 2.6482). The active-surface selection or the combined return of `plasti_rule` differs; found while verifying the per-model visco data (D2), which therefore uses single-model tests.
- **Severity:** LOW. **Effort:** medium (multi-surface return).

### C115. `group_materi_plasti_bounda` wall detection broader than the Professional; missing wall laws and per-index factors (C114 verifier, 2026-11-06)

- **Found:** `group_materi_plasti_boundary_evaluate()` (`group.cc`) marks an element as on a wall when ANY of its nodes is velocity-bounded, whichever `bounda_dof` bounded it, and a listed value also matches element groups (legacy semantics). Professional 02-08-2026 (C114 verifier probe mc_bloc: `group_materi_plasti_bounda 0 3`, `bounda_dof 3 5 -velx`): only element 2 (node 5) is reduced (-25.5416 / -29.9571); the GNU reduces both elements (-20.7494 / -32.471).
- **Also:** the manual lists wall reductions for hardsoil, camclay (M), hypo (deviatoric increments) and the incremental tension (sigy, `_element_group`); the GNU has none. `group_materi_plasti_bounda_factor` takes one value; the manual allows one per listed index.
- **Severity:** LOW. **Effort:** medium (the detection change alters existing output of inputs listing a bounda index whose nodes are not the only bounded ones).

### C114. Wall friction reduction: the GNU does not reduce the cohesion, and the direct Mohr-Coulomb paths ignore or default the factor differently (F5 T3 C2, 2026-11-06) — RESOLVED 2026-11-06

- **Found:** manual 6.687: for `group_materi_plasti_bounda` (and 6.699 `_element_group`, 6.733 `_mpc`) the friction angle phi, the dilatancy angle phi_flow and the cohesion c are reduced by the factor (default 2/3). Professional 02-08-2026 probes (t3/c2, two quad4, Mohr-Coulomb 0.5 1. 0., element 2 on the wall): node 5 sigyy -1.87260 = a group with phi and c both times 2/3 given explicitly (-1.87260); the GNU gives -2.81059 = phi reduced, c kept (-2.80894 in the Professional). `plasti.cc` reduces only phi and phi_flow (five sites: Mohr-Coulomb, its hardening/softening variant, Drucker-Prager, Matsuoka-Nakai and one more).
- **Direct paths:** `materi_direct_full_mc` (no plane normal) applies only the `_wall` records and ignores the factor (Professional: 2/3 by default, `mohr_coul_direct_mpc_total1` input: sigxy 0.0998 -> 0.06662); `materi_direct_cutoff` (plane normal) uses the factor with default 1 instead of 2/3.
- **Severity:** LOW. **Effort:** small code; it changes the existing output of every input with a wall reduction. The manual and the Professional agree (S2), so no owner decision is needed for the cohesion; the defaults of the direct paths also follow the manual.
- **Resolution:** `plasti.cc` reduces c with phi and phi_flow at the five Mohr-Coulomb-type sites; `stress.cc` `direct_wall_factor()` (bounda or mpc factor, default 2/3) on phi, c and phi_flow in both direct modes when no `_wall` record is given (the Professional applies no factor on top of `_wall` values). Probes: direct modes equal the Professional to all printed digits (sigxy 0.0799210 / 0.0800989), Mohr-Coulomb and hardening/softening -1.87283 (Professional -1.87260), Drucker-Prager -2.83937 (-2.83475); Professional `other/mohr_coul_direct8` node_dof identical. Corpus inputs with a wall reduction: only that one and the two `heat_exchanger_pile` timeouts. Tests `group_materi_plasti_bounda_direct1`-`3`, `group_materi_plasti_bounda_druck_prag1`; the `group_materi_plasti_mpc` targets moved to the Professional values. 397a3cb

### C113. `materi_plasti_hypo_substepping -no` has no effect on the Masin routes (F5 T3 C1, 2026-11-06) — DEFERRED (class B)

- **Found:** Professional 02-08-2026, `other/hypo7` (Masin clay): `-no` moves sigyy -231.8090 -> -231.8126. The GNU Masin kernels integrate with a fixed Runge-Kutta substep count (`MASIN_BASIC_NSUB`), so the switch does nothing there (-231.789 either way). The von Wolffersdorff route matches the Professional (-863.7163 both).
- **Severity:** LOW. **Effort:** unknown (the Professional Masin substepping scheme is not documented). DIVERGENCES HYPO-SUBSTEPPING-MASIN.

### C112. `materi_plasti_maximum_iterations` default: manual 20000, GNU 1000 (F5 T3 B, 2026-11-06) — RESOLVED 2026-11-06

- **Found:** manual 6.804: "Default at maximum 20000 iterations are allowed". The GNU cap without the record is `MAX_ITER` 1000 (`stress.cc`), 10 for viscoplastic and 2 for nonlocal viscoplastic groups; `MAX_ITER` also bounds the membrane plane-stress loop.
- **Note:** raising the default only changes runs whose plastic loop does not converge within 1000 iterations (the cap is silent), but those runs change output and time.
- **Severity:** LOW. **Effort:** small code (a separate constant for the plastic loop); owner decision (follow the manual or keep 1000).
- **Decision** (owner 2026-11-06: "vamos a 20000. Es un límite que si hay convergencia más rápida no llega a cumplirse y evitarmos cortes antes de convergencia") and **resolution:** `MAX_PLASTI_ITER` 20000 (`stress.cc`) for the plastic loop; the membrane loop keeps `MAX_ITER` 1000. Impact (instrumented corpus, 931 inputs: test-2026, test-2014, Professional test/other, validation, incremental_driver): only sfnet `examp10` and `examp16` pass 1000 iterations, and both reach 20000 (points that do not converge); every other input is unchanged. `examp10` moves 3e-7 in the displacements; `examp16` node 26 sigyy -218129 -> -237309 (Professional -271692; summed |GNU - Professional| over node_dof 3.61e6 -> 3.44e6), run time 1 s -> 3 s. Test `materi_plasti_maximum_iterations_default1` (examp16 copy, node 26 target).

### C111. `post_point_eps_iso` and `post_point_node_type -node` have no effect in the Professional (F5 T2 V9, 2026-11-04) — DECIDED 2026-11-05 (class A)

- **Found:** Professional 02-08-2026, v9t probes on a 3-element quad4 mesh: a post point 0.0015 beyond the mesh is not found with `post_point_eps_iso 0.1`; with `mesh -follow_material` (edge moved 3 -> 3.1) `post_point_node_type -node` gives the `-node_start_refined` result and a point at 3.05 is not found.
- **GNU:** follows the manual (6.957, 6.959): found with 0.1 (velx 0.1000071), `-node` searches the deformed mesh (tests `post_point_search1`, `post_point_node_type1`). DIVERGENCES POST-POINT-SEARCH.
- **Severity:** LOW. **Effort:** small either way; owner decision (follow the manual or the Professional).
- **Decision** (owner 2026-11-05: "manual+GNU, creo que es lo mejor"): the GNU keeps the manual behaviour; deliberate divergence DIVERGENCES POST-POINT-SEARCH (class A).

### C110. `post_node_rhside_fixed` / `_free` are maxima in the GNU, sums in the manual (F5 T2 V9, 2026-11-04) — RESOLVED 2026-11-05

- **Found:** manual 6.948/6.949: "filled with the sum of node_rhside" for the prescribed / free dofs; the Professional 02-08-2026 gives signed sums (v9pro r0: fixed 2.2e-14, free -2.1e-14). The GNU `post_node_rhside_fixed_free()` (`post.cc`) stores the largest absolute value per dof (2.38 / 1.56e-12 on the same input), so `post_node_rhside_ratio` differs too (v9t t3, von Mises, 2 iterations: Professional 12.94 by default and 1.24 with `-post_node_rhside_free`; GNU 11.49 with the method). The manual 6.950 rule "fixed below 1.e-10 -> the free value" is a threshold 1.e-6 with no fallback in the GNU.
- **Note:** the ratio is read by the iteration control (`top.cc`), so a change moves convergence decisions of existing runs.
- **Severity:** LOW. **Effort:** small code, large blast radius (iteration counts); owner decision.
- **Decision** (owner 2026-11-05: "converjamos. Registramos la suma") and **resolution:** `post_node_rhside_fixed_free()` now stores the signed sums per principal dof. The ratio is the measured Professional rule: per principal dof, the mean |node_rhside| over the free nodes divided by the mean over the fixed nodes (the mean free value when the fixed one is below 1.e-10, and always with `-post_node_rhside_free`), the largest over the dofs. `top.cc` re-evaluates it after the post-solve material pass, so the reported records are those of the final residual (this also resolves C63). Professional probes (c110pro, v9t t3 input): fixed -3.98859 0.0908132 / free 3.98859 -0.0908132, ratio 12.9442 (GNU 12.9430) and 1.23884 with the method (GNU 1.23878). Tests `post_node_rhside_ratio1`, `post_node_rhside_ratio_method1`. The iteration control still uses the per-iteration ratio, with the new formula: `control_timestep_iterations_automatic_default1` 843 -> 723 printed iterations (Professional 654), `_minimum_maximum_wished1` 944 -> 968 (Professional 682); final ratio closer to the Professional in both.

### C109. Generic `post_calcul` item names differ from the Professional `post_calcul_label` (F5 T2 V8, 2026-11-04) — DECIDED 2026-11-04 (class A)

- **Found:** the `.dbs` record `post_calcul_label` (manual 6.904, now written by `print_db.cc`) of the Professional 02-08-2026 names the generic operators with an underscore: `-a_ept -va0_ept -va1_ept -va2_ept -s_ept`. The GNU names have no separator: `-aept -va0ept -va1ept -va2ept -sept`. The pressure, safety and apparent names agree (`-to_pres -st_pres -dy_pres -safety_piping -safety_lifting -young_apparent -poisson_apparent`).
- **Manual:** its example is `-aept` (post_calcul_label names, `target_item` section), so by S2 the GNU keeps its names (DIVERGENCES POST-CALCUL-LABEL-NAMES).
- **Severity:** LOW. **Effort:** small (the name builder in `calcul.cc`, `calculate()`), but renaming changes existing GNU output: vtk / gmsh / node print names and `target_item` labels; owner decision.
- **Decision** (owner 2026-11-04: "mantenemos nombres GNU-manual"): the GNU keeps the manual names; deliberate divergence DIVERGENCES POST-CALCUL-LABEL-NAMES (class A).

### C108. V7 per-model dof labels and total strain kappa measures differ from the Professional (F5 T2 C107 probes, 2026-11-04) — RESOLVED 2026-11-04

- **Correction of V7:** the V7 claim "the Professional accepts the initialisations but adds no dof" is wrong. The probes of `0fa067b` (scratch `v7pro`, kappa_shear_ip_state1 model, von Mises) show that the Professional 02-08-2026 adds every dof except `_camclay`, which it rejects.
  - **Per-model plastic strains:** values equal to GNU (eppvo* = epp* under von Mises, the other models 0; C107 probe on druckprag tria3: identical to 6 digits). Labels differ: Professional `eppvo` / `eppmo` / `eppte` / `eppca` / `eppco` / `eppdi` / `eppdr` / `eppha`, GNU `eppvm` / `eppmc` / `epptn` / `eppcap` / `eppcmp` / `eppdip` / `eppdrp` / `epphs`. The GNU cap yy label is the typo `eppcappyy`.
  - **Total strain kappas:** labels `kapto` / `strtocokap` / `strtoshkap` / `strtotekap` (GNU `eptk` / `eptkc` / `eptks` / `eptkt`) and different values. Node 1: kappa 0.0076321 vs GNU 0.008065; shear 0.007231 vs 0.0080105; compression -0.0100122 (signed) vs GNU +0.0100122; tension 0 vs GNU 0.0040283.
- **Manual:** 4.45-4.48 give no formula ("maximum strain size", "maximum principal compression total strain"), so by S2 the Professional is the oracle for the measures.
- **Docs to rewrite with the fix:** DIVERGENCES STRAIN-PLASTI-MODEL-DOFS, the SEGUIMIENTO V7 row and lines, the V7 user and developer pages.
- **Severity:** LOW. **Effort:** medium (measures to identify; renaming the older labels changes existing GNU output, owner decision).
- **Labels resolved 2026-11-04** (owner: "si, pásalas"): `database.cc` uses the Professional labels for the 8 per-model plastic strains (`eppca`, `eppco`, `eppdi`, `eppdr`, `eppha`, `eppmo`, `eppte`, `eppvo`; the `eppcappyy` typo is gone) and the 4 kappas (`kapto`, `strtocokap`, `strtoshkap`, `strtotekap`). Camclay keeps `eppcc` (the Professional rejects that initialisation). The `mstrain_*` targets were renamed; RED: the old binary rejects `-eppdrxy` ("I do not know").
- **Measures resolved 2026-11-04** (S2): analytic Professional probe on `materi_strain_total_kappa1`. The kappa equals GNU. The compression kappa is now signed (negative, history minimum, `materi.cc`). The Professional's current (non-history) compression / shear values and its zero tension kappa contradict the manual ("as occurred in history"): DIVERGENCES STRAIN-KAPPA-HISTORY. The V7 docs (DIVERGENCES, SEGUIMIENTO, manual pages) were corrected. The remaining non-uniform differences (node 1 kapto 0.0076321 vs 0.008065) come from the total strain itself (ept differs in that von Mises run), not from the kappa measure.

### C107. The older per-model plastic strain dofs ignore the model and the epp recovery (F5 T2 V7, 2026-11-03) — RESOLVED 2026-11-04

- **Found:** `materi_strain_plasti_hardsoil / _cap / _compression / _diprisco / _druckprag` (`ddacc43`, `materi.cc` consolidated block) get the plastic strain increment in every element group, whatever its model.
  - **The Professional:** 02-08-2026, kappa_shear_ip_state1 model (von Mises): `_druckprag`, `_hardsoil` and `_cap` dofs are added and stay 0.
  - **The GNU:** these dofs use the h_lump weight, and `general.cc` lists none of their types in the strain recovery weight list. On tria3 their nodal value differs from epp at shared nodes (node 2: eppdrp -0.00164 vs eppxx -0.00166).
  - **V7 dofs:** the 4 new dofs (camclay / mohr_coul / tension / vonmises) are filled only for the group's model, with the epp recovery.
- **Severity:** LOW. **Effort:** small. Move the 5 dofs to the V7 block. This changes existing output (the `mstrain_*` tests), so it needs an owner decision first.
- **Resolution** (owner 2026-11-04: "sí, implementalo así. Los llevamos a cero si no tienen ese modelo"): the 5 dofs moved to the per-model block of `materi.cc` (model gating with `get_group_data`: hardsoil, cap / cap1, compression / compression_direct, diprisco, druckprag) and to the `general.cc` strain recovery weight list; the old h_lump block was removed. Under another model they stay 0, in their model they equal epp. RED: `materi_strain_plasti_model5` (druckprag on tria3) gave eppcap* = epphs* = 0.913406 with the old binary. The `mstrain_*` and `mhardsoil_*` targets still hold.

### C98. Meaning of `print_vtk_coord` / `control_print_vtk_coord` (F5 W2 V2, 2026-11-02) — RESOLVED 2026-11-02

- **Measured** (Professional 02-08-2026): every vtk file has `POINTS` and a
  point field `VECTORS coord`; `print_vtk_coord -no` or
  `control_print_vtk_coord 0 -no` removes only that field, a control record
  overrules the global one. The GNU has no `coord` field and its `-no`
  (Sprint 11 decision, manual 6.340 "the coordinates of nodes is plotted in
  vtk") omits `POINTS`. The F5 global record keeps the GNU meaning.
  DIVERGENCES `C98` (class B).
- **Severity**: LOW. **Effort**: small (add the field, change the `-no`
  branch; owner decision first: it adds a field to every GNU vtk file and
  changes an existing record).
- **Resolution** (owner 2026-11-02 "dale, seguí a Pro"): `print_vt.cc`
  always writes `POINTS` and adds `VECTORS coord` (stored node coordinates)
  first in `POINT_DATA`; `-no` drops only that field. Tests test-2014
  `vtk_coord1`, test-2026 `print_vtk_coord1/2` updated. DIVERGENCES row
  removed.

### C97. Professional `tutorial_4` does not run in the GNU; its nonlinear result differs (F5 W1 U6, 2026-11-02) — OPEN

- **Measured** (`test/tutorial/tutorial_4/tutorial_4.dat` with
  `test_calculation true`, the `linear_calculation_apply` example): the GNU
  stops with "post_calcul -materi_stress -force is a NODAL calculation;
  POST_LINE_DOF/POST_POINT_DOF/POST_QUADRILATERAL_DOF records are not
  supported" (the input has `post_point`), and without `-force` with
  "Error detected in gid 6 printing" (`control_print_gid -separate_index`).
  Without both: `post_node_result 0` velx -2.29 GNU, -4.17 Professional
  02-08-2026 (`linear_calculation_apply -no`); with `-yes` the GNU ends at
  `post_node_rhside_ratio` 4.6e-13, the Professional at 3.9e-2.
- **Severity**: LOW. **Effort**: medium (three separate causes: the
  `-force` + `post_point` refusal, the GiD print, the excavation /
  interface / groundflow result).

### C96. Integration-point state of an element added by `control_input` (F5 W1 U6, 2026-11-02) — OPEN

- **Measured** (`control_input2`: one quad4, a second quad4 and two nodes
  added at index 1, two steps): the nodes of the new element agree (sigyy
  20), the shared nodes average to -10 in the GNU and -5 in the
  Professional 02-08-2026. The GNU starts the new element from the node
  values (`element_dof_initialised` 0, as for the elements of the mesh
  operations); the Professional rule is not identified. DIVERGENCES
  `CTRL-INPUT-NEWELEM` (class C).
- **Severity**: LOW. **Effort**: small/medium (measure the Professional
  rule for new elements, then decide).

### C95. `mrange` and the expanded length of value ranges (F5 W1 U6, 2026-11-02) — OPEN

- **Measured**: the Professional 02-08-2026 checks `mrange` against the
  words of a value range (`bounda_dof` node list) and against the words
  plus markers and the expanded list of an index range; the GNU does the
  same at input. The expanded length of a value range (`-from a -to b` in
  a node list, expanded at run time in `bounda.cc`, `data.cc`,
  `filter.cc`) is checked in the GNU only against the compiled `MRANGE`;
  whether the Professional checks it against `mrange` was not measured.
- **Independent review (2026-11-02)**: the first U6 probe notes record
  that the Professional does not limit a range inside the values of a
  record such as `geometry_set`, while the GNU checks every integer value
  range of the data part against `mrange` (input.cc); only relevant when
  `mrange` is given. `mrange` / `mstring` written in the data part are
  still rejected, now with "Unexpected end of input" instead of "I don't
  know what to do with" (the names are registered for the initialization
  part). To measure together with the expanded length.
- **Severity**: LOW. **Effort**: small (one probe, then a check in the
  three run-time expansions if needed).

### C94. Automatic time stepping: step-size law differs from the Professional (F5 W1 U4, 2026-11-02) — DEFERRED (class B)

- **Measured** (4x4 quad4 von Mises, `control_timestep 0 0.05 1.`,
  `control_timestep_iterations_automatic 0 1.e-4 1.e-6 0.5`, iterations
  from `control_print_number_iterations`): both codes now stop a step
  after 8 iterations by default (manual 6.387, 2 8 4). The Professional
  grows the step by 1.25 (sometimes 1.5625) after a step that converged
  in fewer iterations than wished and divides it by 4 after a step at the
  maximum; the GNU multiplies by 1 + (wished - iterations)/20 and divides
  by 10. 166 Professional steps against 239 GNU steps to t = 1.
  `control_timestep_adjust_minimum_iterations -yes` (Professional default)
  has no specified rule; the GNU keeps the minimum.
- **Corpus effect of the 2 8 4 default** (omen `u4_runs`, the 17 inputs
  with automatic stepping, before / after `9d0da65`): 14 identical; sfnet
  `truss7` now reaches the maximum, the step is divided by 10 until it is
  below the minimum and the index ends silently at t = 0.069 instead of 1
  (static load: the target disy -1 is still met, rc 0); sfnet `taylor3`
  differs (node_dof up to 1.3 relative; its target fails before and after);
  Professional other `vajont_slope_fs` stopped before at the minimal time
  step (t = 0.16, rc 1) and now runs past 900 s. The Professional has no
  reference for the three (parse error or timeout). A step-size collapse
  ends the index without a message unless
  `control_timestep_iterations_automatic_stop -yes`.
- **Severity**: LOW (results converge; the step sequence and cost
  differ). **Effort**: medium (measure the Professional law on several
  iteration patterns).

### C93. `post_element_force` is 0 with `global_element_dof_apply -no` (F5 W1 U2, 2026-11-02) — OPEN

- **Measured**: quad4 of `control_post1` pulled in y, 5 steps,
  `global_element_dof_apply -no` (= `options_element_dof -no`): GNU
  `post_element_force_result 0` = 0, Professional 02-08-2026 3.2967033
  (the same as with `-yes`; the nodal stresses agree). The GNU evaluates
  the section forces from `ELEMENT_DOF`, which `-no` does not store.
- **Severity**: LOW. **Effort**: small/medium (fall back to the nodal
  stresses, or reject the combination).

### C92. Thread overhead; `processors_maximum` default deferred (F5 W1 U2, 2026-11-02) — DEFERRED (class B, F6)

- **Measured**: `parallel_sys_routine()` creates and joins `nthread`
  pthreads per call and `element_loop()` allocates `nthread` copies of
  the nodal arrays per sweep. `validation1_steady1` (80 elements, 300
  steps): 12.0 s on 1 thread, 15.5 s on 4, 65.5 s on 16 (83 s system
  time; it hits the build_safe timeout). A 22 500-element plastic quad4
  model: 23.2 / 13.7 / 13.5 s (the serial solver dominates beyond 4).
- **Decision (owner 2026-11-02)**: `processors_maximum` is implemented
  but the GNU default stays `-no` (Professional `-yes`) until the thread
  pool of F6; DIVERGENCES `PROC-MAX`.
- **Found with it and fixed**: `force_point` was added once per thread
  (`force_point_threads1`).
- **Corpus A/B, same binary, 1 vs 4 threads (omen `mt_runs`, after the
  force_point fix)**: tests whose final `.dbs` differ by more than 1e-6
  relative: Professional other 60/330 (`interface_patch` and `patch1`
  abort with signal 6; node / post values off by > 1 % in `force2`,
  `incnav5`, `ground11_phreatic_level`, `distri3`, `ground6`, `force4`,
  `ground4`, `excavate1`, `ground8`, `mohr_coul_direct4`, `bimet1`,
  `matrix3`, `ground11_nonsaturated`), sfnet 41/197 (`taylor3`,
  `examp16`, `force2`, `examp17`, `incnav5`, `examp18`), test-2026 50/294
  (the `reduced_*` family and `spring_oblique_free1` in the residual
  norms), test-2014 6/244; the rest are round-off or residual-norm
  (`post_node_rhside_*`) changes. Multithreaded runs are therefore not
  reliable yet: F6 must make them equal to the single-thread results
  before `processors_maximum` can default to `-yes`.
- **Independent review of `1a4cb7d` (pre-existing thread hazards, for
  F6)**: `parallel_contact` calls `area_node_dataitem()` (contact.cc:141)
  once per thread, an idempotent but redundant and unlocked whole-mesh
  sweep of database writes (it can abort with "Data is allocated in a
  parallel loop" when a target record is new); `parallel_new_dof_before`
  PUTs `NODE_INERTIA` (VERSION_NEW) inside the parallel region;
  `parallel_post_node` (post.cc:294-330) re-partitions through the shared
  `inext` counter in a `for(;;)` loop (fragile coverage with more than one
  thread). A restart from a `.dbs` carries the `options_processors` value
  of the machine that wrote it.
- **Severity**: LOW. **Effort**: medium (F6: persistent pool, buffers
  allocated once, then flip the default).

### C90. `control_print_history` file names differ from the Professional (F5 W1 U3, 2026-11-02) — OPEN

- **Measured**: `control_print_history 0 -post_point_dof 0 -disy` writes
  `disy0.his` in the GNU and `post_point_dof_0_disy.his` in the
  Professional 02-08-2026; `-post_element_force_result 0 0` writes
  `00.his` / `post_element_force_result_0_0.his`. The contents (time,
  value) are equal; both append to an existing file.
- **Impact**: scripts reading the Professional history files do not find
  the GNU ones; two GNU history items with the same dof and index (e.g.
  `-post_point_dof 0 -disy` and `-node_dof 0 -disy`) may share one file.
- **Severity**: LOW. **Effort**: small (file name builder of the history
  print; check the corpus scripts that read `.his` names).

### C91. A switch value read into a DOUBLE PRECISION record is stored as a raw enum (F5 W1 U3, 2026-11-02) — OPEN

- **Measured**: `force_edge_water 10 -yes` (test-2026
  `phreatic_static_submerged1`) is accepted into the DOUBLE PRECISION
  record `force_element_edge_water` and printed in the `.dbs` as
  `force_element_edge_water 10 -1318`, the enum number of `-yes`; the
  number changes whenever an enum is added (-1312 before F5 W1 U1/U3).
  The calculation is not affected, but the record cannot be read back
  (restart / `.dbs` reuse) and the manual record syntax should be checked.
- **Severity**: LOW. **Effort**: small (check manual 6.5xx syntax of
  `force_edge_water`, reject or convert the switch).

### C87. Free elastic interface in traction: default gap (F4 T6f C60, 2026-10-19) — RESOLVED 2026-10-21 (owner decision: follow the manual and the Professional)

- **Fix (F4 T6f C87, 2026-10-21)**: `interface_element()`
  (`interface.cc`): default gap 0 instead of 1.e20; an opened interface
  (gap or tension limit) skips the Mohr-Coulomb return, has zero shear
  force (record, rhs, history) and `kt * residual_factor` in the
  matrix. Probes with blocks attached (C60 geometry, four node orders,
  kn = kt = 1000, pulled 0.01, slid 0.002): no gap record Professional
  0 / 0 (GNU before 10 / +-2), gap 0 Professional 0 / 0 (GNU before 0 /
  +-2), gap 1.e20 both 10 / +-2, compression both -10 / +-2. The manual
  agrees (6.625 "If you want to allow always tension stresses in an
  interface set gap to, by example, 1.e20"; 6.632 an opened interface
  "does not have stresses"). Tests `interface_tension_default_gap1`,
  `interface_open_gap0_shear1`, `interface_tension_gap1e20_1`;
  `iface_mc_dil_1step` target -1.49338 -> -2.97197 (= Professional),
  `iface_mc_dil` gets an explicit gap 1.e20. Found C89.


- **Measured** (Professional 02-08-2026): one quad4 interface, no
  continuum element, kn = kt = 1000, no `group_interface_gap`, side 2
  moved by 0.001 normal to the interface for the node orders `1 2 3 4`,
  `3 4 1 2`, `2 1 4 3`, `4 3 2 1`. Professional normal stress -1 / 0 /
  0 / -1, GNU -1 / +1 / +1 / -1: the orientation agrees (C60), but the
  Professional carries no tension on the opened interface.
- **Cause (hypothesis)**: the GNU default gap is `1.e20` (always
  closed, `interface_element()`); the manual sentence "if you want to
  allow always tension stresses in an interface set gap to, by example,
  1.e20" suggests a Professional default near 0. Not yet checked with
  blocks attached (the C60 block probes were all in compression).
- **Decision needed**: follow the Professional (default gap 0) changes
  every interface in traction without a gap record. **Severity**: LOW.
  **Effort**: cheap/medium (plus corpus A/B).

### C86. `post_element_force` with truss / beam / spring / interface elements in the selection (F4 T6f LOW queue Q5, 2026-10-17) — RESOLVED 2026-10-19 for trusses and springs (owner decision: follow the Professional)

- **Fix (F4 T6f C86, 2026-10-19)**: `post_element_force_calculate()`
  (`calcul_force.cc`) adds the nodal force of `-truss`, `-spring1` and
  `-spring2` elements (the `element_rhside` of `truss.cc` / `spring.cc`,
  from `ELEMENT_TRUSS_FORCE` / `ELEMENT_SPRING_FORCE`) on the side of
  their centroid. Professional 02-08-2026 probes, all equal: truss on a
  boundary section 100 / 10 / -40, two trusses at an interior section
  node 100 / 10 / -20 (continuum alone 55 / 10 / 2.5), spring2 along the
  section normal 100 / 10 / -40 (spring 55.385), diagonal truss 1-10 100
  / 10 / -40, 3D hex8 strip + truss 100 / 10 / -4 / -40 / 16. Tests
  `post_element_force_2d_truss1`, `_truss_interior1`, `_spring2_1`,
  `post_element_force_3d_truss1`; `post_element_force_2d_beam_error1`
  keeps the beam refusal. A spring1 on a section node counts on the
  positive side (interior probe 100 / 5.496711 / -20); the spring1 and
  oblique-spring section tests were added after C88 fixed the spring
  forces (`post_element_force_2d_spring1_interior1`,
  `post_element_force_3d_truss_spring1`). Not followed: beams and truss-beams (the
  Professional 2D beam carries no force in the probe, no reference;
  BACKLOG C37), contact springs and interfaces (not measured): they stay
  refused, DIVERGENCES C84-REFUSE.

- **Measured** (Professional 02-08-2026, `scratchpad/lowq/q5/`): the C84
  quad4 strip (section x=0, loads 50 / 5 at x=4) plus a `-truss` element
  (E A 1e4) from node 1 (on the section) to node 5: the Professional
  INCLUDES the truss in the section (normal 100 = truss 45 + continuum
  55, shear 10, moment -40 = 45 x (0 - 0.5) + continuum -17.5): the
  analytically correct total section force. A truss that does not touch
  the section (node 2 -> 5) is not counted (100 / 10 / -40 from the
  continuum alone). With `post_element_force_group` excluding the truss
  both codes give 55 / 10 / -17.5. The GNU stops (`Error:
  post_element_force does not support element : 5`).
- **Not measured**: beams, springs (`-bar2` springs are refused by the
  Professional in 2D), contact springs, interfaces.
- **Decision needed**: the Professional neither ignores nor mistreats the
  truss, so the plan rule (follow if ignored, stop if wrong) does not
  apply: implementing the inclusion means adding the structural-element
  nodal forces (truss: N along the axis; beams: shear and moments) of the
  elements with a node on the section. DIVERGENCES C84-REFUSE.
- **Severity**: LOW. **Effort**: medium (truss cheap, beams/interfaces
  need probes).

### C84. `post_element_force` 2D/3D section on non-tensor elements and in 2D (F4 T6f C66, 2026-10-16) — RESOLVED 2026-10-17 (`f2ae525`)

- **Follow-ups (F4 T6f LOW queue, 2026-10-17)**: the 1D result is one
  value, as the Professional `force16` (`feafe50`); the refusal of
  non-continuum elements is now item C86.

- **Symptom**: the 2D/3D `post_element_force` path of `calcul_force.cc`
  (the element scan near line 2800) chooses `npol` 2 for
  `quad4/tria3/hex8/bar2` and 3 otherwise, builds the rule with
  `msf_element_rule()` (a tensor-product rule) and integrates the
  ELEMENT_DOF stresses with the quad/hex helpers
  `msf_element_internal_forces_2d/3d()`, whatever the element type. A
  tria3, tria6, tet4, tet10, prism6 or prism15 therefore uses
  quad/hex shape functions, a 2x2(x2) / 3x3(x3) point count that is not
  the element's, and reads ELEMENT_DOF slots past its real points.
  Probe (2026-10-16, `scratchpad/c66/pef/`): `tpef_force10` (cantilever,
  load (0,1000,-123.4), section at the fixed face) with each hex8 split
  into 6 tet4: GNU normal 3143.4 / shear1 -539.4 / moment1 -16896;
  Professional 02-08-2026 -123.4 / 1000 / -100000 (= the hex8 of both
  codes). 2D: the GNU requires 7 values in the record while the
  Professional accepts 6 (`n_x n_y s0_x s0_y mid_x mid_y`) and rejects
  a 7th; the GNU skips the 6-value record ("invalid length") and with
  7 values prints NaN for every result, quad4 included. The
  Professional gives normal 100 / shear 10 / moment -40 for a 4x1
  cantilever strip of quad4 and of tria3 (exact by equilibrium).
- **Cause**: the section device was written and measured for
  quad4/quad9/hex8/hex27 (the support L5 helpers) and its 2D record
  length was never checked against the Professional.
- **Next step**: per element type, integrate B^T sigma with the
  element's own shape functions and `pol()` rule (or reuse the element
  internal-force assembly); accept the 6-value 2D record; find the 2D
  NaN source; tests from the probes above.
- **Severity**: LOW (post-processing only, wrong numbers without a
  warning). **Effort**: medium.
- **Resolution** (F4 T6f C84, `f2ae525`): manual 6.927 gives the 2D
  record as 6 values (no input of the corpus, test-2014, test-2026 or
  sfnet used 7; the Professional `force14`/`force15` use 6), so the
  GNU follows it: `input.cc` stops in the data part on a wrong count
  (`Number of data values expected : 6`, 12 in 3D, 2 in 1D, as the
  Professional) and the 2D result has the Professional 3 values
  (normal, shear0, moment0). The 2D NaN was the stress copy (`sxx sxy
  sxz` in a 6-stride array read as `sxx syy sxy` with stride 3). The
  non-tensor continuum elements (tria3, tria6, quad16, tet4, tet10,
  prism6, prism15) and every 2D element of an axisymmetric group
  integrate with `pol()` and the `ELEMENT_DOF` stresses of their real
  points (`pef_element_internal_forces_pol`; axisymmetric: ring forces
  with the hoop term); quad4/quad9 plane and hex8/hex27 keep the
  tensor helpers (hex8/hex27 byte-identical). With elements on both
  sides of the section (no `post_element_force_normal`) the result is
  the mean of the positive side and the negated negative side (the
  plain sum cancelled; Professional corpus `force14`/`force15` now
  equal). Other elements are
  refused with an error (Professional not measured: DIVERGENCES class
  B C84-REFUSE). Professional = GNU on every probe (tet4, tet10,
  prism6, prism15 cantilever -123.4 / 1000 / -100000; quad4, tria3,
  quad9, tria6 strip 100 / 10 / -40; axisymmetric quad4, quad9, tria3,
  tria6 ring strip, normal 100; interior strip section 100 / 10 / -20;
  reversed normal). Tests `test-2026/post_element_force_*` (14, RED on
  `54b8170`).

### C83. tria3 / tet10 `-maximal` point counts and the axisymmetric tria3 `-maximal` (F4 T6f C67 remainder, 2026-10-16) — RESOLVED 2026-10-16 (part 2, `72f78bd`; part 1 decided: deliberate divergence from the manual, class A)

- **Symptom**: manual 6.623 specifies "a four-point integration scheme"
  for the tria3 `-maximal` and "a ten-point scheme" for the tet10
  `-maximal`; the Professional 02-08-2026 uses 3 points at the nodes and
  5 points (Keast), and so does the GNU (no difference between the codes
  in a plane group). In an axisymmetric group the GNU `pol()` forces the
  1-point rule for every tria3; the Professional keeps the 3 nodal points
  for `-maximal` and switches to 1 point only in an element with a node
  on the axis (like C71); nodal stresses of the probe 3 % of the scale
  apart.
- **Decision** (owner 2026-10-16): part 1, keep 3 / 5 points (GNU =
  Professional), deliberate divergence from the MANUAL: the manual gives
  no weights and the classic 4-point degree-3 triangle rule has the
  negative -27/48 centroid weight rejected in C82. A positive 10-point
  degree-3 tet rule exists (Keast 1986: (a,b,b,b) a=0.5684305842 weight
  0.2177650699 + 6 edge midpoints weight 0.0214899534, normalised to 1)
  but no version uses it; the 10-node nodal rule is not positive (vertex
  weight -1/20). Part 2, follow the Professional.
- **Resolution** (part 2, `72f78bd`): `pol()` no longer forces 1 point
  for an axisymmetric tria3; `-maximal` keeps the 3 nodal points and an
  element with a node at |r| < EPS_COORD gets the 1-point centroid rule
  (Professional probes: the same per-element node-radius switch, as C71 /
  C80). Why the GNU forced 1 point: the line is in the 2016 initial
  import without a comment; the only reason found is the r=0 hoop
  singularity (h/r infinite, volume 0 -> NaN, C71), which the on-axis
  switch avoids. Default and `-minimal`: 1 point in both codes
  (unchanged). Recovery unchanged (tria3 1/3-point branch of
  `sri_stress_recovery_weight()`). 27 tria3 probes (plane, axisymmetric
  off/on the axis, 3 methods x 3 point values) = Professional within
  3.1e-12 (before 3.2e-2). Test `test-2026/recovery_tria3_axisym_maximal1`.
- **Severity**: LOW (non-default options). **Effort**: cheap.

### C82. `prism15` `group_integration_points -minimal` (F4 T6f C67 remainder, 2026-10-16) — DECIDED 2026-10-16 (deliberate divergence, class A)

- **Symptom**: the Professional 02-08-2026 integrates a prism15 with
  `-minimal` with 8 points: the 4-point degree-3 triangle rule (centroid
  weight -27/48, three points 25/48) x 2 Gauss levels, recovered with the
  physical linear fit. The GNU keeps 21 points (the record is not
  branched on); nodal stresses of the probe differ 0.9 % of the scale.
- **Analytic check**: the 8-point element stiffness of a distorted
  prism15 has 3 NEGATIVE eigenvalues (the negative centroid weight); the
  rule is not positive definite, unlike the 21-point rule (6 zero
  eigenvalues, the rigid modes).
- **Decision** (owner 2026-10-16): keep 21 points, deliberate divergence
  (class A; DIVERGENCES C82, SEGUIMIENTO Descarte); docs only.
- **Severity**: LOW (non-default option). **Effort**: cheap (docs) /
  cheap (code, if followed).

### C66. tria6 integration rule differs from the manual and the Professional (F4 T6f, 2026-10-11; re-scoped 2026-10-16) — RESOLVED 2026-10-16 (`175faf1`; lumping: deliberate divergence TRIA6-HRZ, class A)

- **Symptom**: two-tria6 probe (`test-2026/tria6_nodal_recovery_mean1`):
  displacements, integration-point and single-element nodal stresses
  differ 1e-6..2e-5 relative from the Professional 02-08-2026 (the
  nodal recovery itself is the same mean). The tria3, tet4 and bar2
  probes agree to 1e-13.
- **Cause (hypothesis)**: the tria6 6-point nodal rule weights
  (`polynom.cc`, corner weight `EPS_TRIANGLE` 1e-6, mid-side 1/3-1e-6)
  or a different Professional tria6 rule.
- **Cause (measured 2026-10-16, F4 T6f C67 remainder)**: a different
  rule. Manual 6.623: "For -tria6 elements a seven-point scheme will be
  used for -maximal and a four-point scheme will be used with -minimal";
  the Professional does exactly that (default = `-maximal`), recovers
  the points with the physical linear fit (residual 5e-13) and the
  arithmetic mean, and accepts axisymmetric groups (7 / 4 points). The
  GNU uses its 6 nodal points for every value (exact for the degree-2
  stiffness of a straight-edged element, hence only ~1e-5 apart there)
  and refuses axisymmetric tria6. On a curved, distorted 12-element mesh
  the nodal stresses differ 7.7 % (`-maximal`) / 9.5 % (`-minimal`) of
  the scale.
- **Professional lumping (measured 2026-10-16, F4 T6f C66)**: one
  element, density 1, E 1e-9, a unit point force at one node with every
  other node fixed, one step dt 0.01, m = F dt / v (the same method as
  `tet10_mass_hrz1`); `condif` capacity cross-checked with a point heat
  source. The Professional uses the ROW SUM of the consistent mass for
  every element, independent of `group_integration_points`:
  - tria3 A/3, quad4 A/4, quad9 A/36 / A/9 / 4A/9 (corner / mid-side /
    centre), tet10 corner -V/20 and mid-edge V/5 (as recorded in C31):
    the same as the GNU, except the tet10 (GNU HRZ, TET10-HRZ).
  - tria6 plane, default / `-maximal` / `-minimal`: corner 0 (< 1e-11;
    a free single-element dynamic step fails with a Pardiso zero pivot),
    mid-side A/3. Total A.
  - tria6 axisymmetric (r = 1..2, velocity and capacity alike): corner
    2 pi A (r_i/30 - (r_j + r_k)/60), i.e. -pi/60 at the two r=1
    corners and +pi/30 at the r=2 corner; mid-side 1.466 / 1.257 /
    1.466. Total 2 pi A r_c. Negative corner masses away from the axis.
  - GNU HEAD (6 nodal points): corner 5e-7 (the 1e-6 corner weight),
    mid-side A/3; axisymmetric refused.
  Candidates for the plane unit tria6 (A = 1/2): row sum 0 / 1/6
  (Professional), HRZ A/19 / 16A/57 (0.0263 / 0.1404), nodal quadrature
  ~0 / 1/6 (GNU), equal split A/6.
- **GNU lumping sites with an interior tria6 rule**: `general.cc`
  inertia block, `weight = h[inol]` (row sum): the principal dofs
  (`materi_velocity`, `condif_temperature`, `groundflow_pressure`,
  `wave_fscalar`) and the h-weighted state dofs (density, history
  variables, hypo/camclay histories, `plasti_f`, softening variables,
  strain energy, intergranular and model-specific plastic strains, void
  fraction, work) get zero (plane) or negative (axisymmetric) corner
  entries. The recovered dofs (stress via `sri_stress_recovery_weight()`,
  strains, kappa, kappa_shear, cap1, damage, rho, Maxwell via
  `strain_recovery_weight()`) use the per-element fit F/volume and are
  not affected. `calcul_force.cc` lumped shares (1D only),
  `elem.cc` relaxation and `interface.cc` storage are not affected.
- **Gate (2026-10-16)**: stopped before any source edit. Following the
  Professional lumping gives zero (plane) or negative (axisymmetric)
  corner masses, which the gate excludes.
- **Next step (owner decision)**: adopt the manual 7 / 4-point rule with
  the linear-fit recovery (the recovered dofs are unaffected), and lump
  the h-weighted dofs of the tria6 with HRZ (corner A/19, mid-side
  16A/57 plane; r-weighted diagonal scaled to 2 pi A r_c axisymmetric),
  a deliberate divergence like TET10-HRZ; or keep the GNU 6-point rule
  (class A). The tria6 corpus cases must then be re-checked against the
  Professional.
- **Status**: **RESOLVED 2026-10-16 (F4 T6f C66 phase 2, `175faf1`)**,
  owner decision 2026-10-16: the manual rule (7 points default /
  `-maximal`, 4 points `-minimal`) with the Professional physical
  linear-fit recovery and arithmetic mean; axisymmetric tria6 accepted
  (Professional: no rule switch and no h-average in an element with a
  node on the axis, the points are interior, no NaN); HRZ lumping of the
  h-weighted dofs (plane corner A/19, mid-side 16A/57; axisymmetric
  r-weighted diagonal scaled to 2 pi A r_c; the same for every rule),
  deliberate divergence TRIA6-HRZ (class A). Every probe node = the
  Professional within 3.1e-10 of the scale. Tests `recovery_tria6_curved1`,
  `recovery_tria6_minimal1`, `recovery_tria6_axisym1`, `tria6_mass_hrz1`,
  `tria6_mass_hrz_axisym1`. Corpus: omen A/B c66_runs vs c83_runs (997 runs): changed `tria6` (pro_other and sfnet copies; control_mesh_split quad9 -> tria6, uniform stress, node_dof change 3.3e-10, vs Professional 4.67e-5 before and after) and `tria6_nodal_recovery_mean1` (5.1e-6, now = Professional within 5e-10); only-new the 5 tests; `slope_classical_numerical` (584 tria6, TIMEOUT 120 s on both) rerun with 1800 s: collapse factor old 1.26297 / new 1.26265 / Professional 1.27045 (-0.59 % -> -0.61 %, step-control level), final nodal stresses vs the Professional 20.8 -> 0.25 of the column scale (toward the Professional). Follow-up C85 (groundflow
  velocity numerator, condif absorption).
- **Severity**: LOW. **Effort**: cheap/medium.

### C85. tria6 h-weighted lumped terms outside the HRZ conversion (F4 T6f C66 phase 2, 2026-10-16) — RESOLVED 2026-10-17 (F4 T6f C85)

- **Symptom**: two h-weighted lumped terms were not converted to the
  tria6 HRZ weight because their files were outside the C66 work unit:
  (1) `groundfl.cc` groundflow velocity numerator `volume * h * (...)`;
  `general.cc` therefore keeps `h` for `groundflow_velocity`, whose
  tria6 corner denominator `sum volume h` is 0 with the 7 / 4-point
  rules (a tria6 groundflow velocity is singular at the corners; no corpus
  case uses it); (2) `condif.cc` lumped absorption
  (`volume * h * condif_absorption` on the diagonal): ~0 at the tria6
  corners, as with the old nodal rule.
- **Resolution** (owner decision 2026-10-17, F4 T6f C85; class A
  TRIA6-HRZ extended): `groundfl.cc` and `condif.cc` use
  `lumped_hrz_active ? lumped_hrz_weight[inol] : h[inol]` (the
  `materi.cc` `h_lump` pattern) and `general.cc` drops the
  `groundflow_velocity` exception. Tests (hand targets, RED on
  `c6fa9aa`, GREEN): `test-2026/tria6_groundflow_velocity_hrz1` (linear
  pressure, k 10: gvelx 10 at every node; before corners 35.2 / -40 /
  -144 / -368, 0/0 round-off), `tria6_condif_absorption_hrz1` /
  `_axisym1` (dT/dt = -T, one implicit step dt 0.1: 1/1.1 at every node;
  before plane corners 0.9955 / 1.0027 and mid-sides 0.8946 / 0.8934,
  axisymmetric corners up to 1.137). Other elements: `lumped_hrz_active`
  0, unchanged (quad4 control in the absorption tests; corpus A/B).
- **Severity**: LOW. **Effort**: cheap.

### C65. Professional `post_point_dof` at a node differs from its own `node_dof` (F4 T6f, 2026-10-11) — DEFERRED (class B, F4 T6f LOW queue, owner plan 2026-10-17)

- **Deferred (2026-10-17)**: Professional-side post-processing inconsistency, folded into
  DIVERGENCES PPD (the Professional evaluates post points from the
  element-local recovered value); no GNU change planned.

- **Symptom**: Professional `validation_14`: `post_point_dof` at the
  coordinates of a node gives `eppxy` 0.00815 while its `node_dof` at that
  node is 0.01400. GNU interpolates the nodal field, so its post point at
  a node equals its `node_dof`. Consistent with the known AUDIT section 5
  item (the Professional evaluates post points from the element-local
  recovered value, `element_dof`); which element value (or combination)
  it uses at a shared node is not identified.
- **Next step**: identify which element the Professional uses (lowest
  number, first found) on a two-element probe; decide whether to follow.
- **Severity**: LOW (post-processing only). **Effort**: cheap.

### C64. quad4 per-element recovery of the Professional: linear `1,xi,eta` least squares? (F4 T6f, 2026-10-11) — RESOLVED 2026-10-12 (F4 T6f K21 WU1)

- **Resolution**: the Professional fits `1,x,y` in PHYSICAL coordinates
  (quad4 2x2, quad9 2x2; 1-point quad4 = point value; axisymmetric
  element with a node at r=0: `sum h sig / sum h`), measured exact on
  `tube1`, `examp1`, `incnav4` and distorted probes. Implemented per
  element in `miscel.cc` `recovery_fit_factor()`; tests
  `test-2026/recovery_*`. The natural `1,xi,eta` fit matched 92/95
  `validation_14` elements because it equals the physical fit on
  (near-)parallelograms. Remaining types: C70 (resolved 2026-10-13).

- **Symptom**: on a Drucker-Prager variant of `validation_14` (quad4,
  2x2 Gauss) the Professional per-element nodal values (`element_dof`
  from `element_intpnt_dof`) are reproduced exactly in 92 of 95 elements
  by a linear least-squares fit `1,xi,eta` of the 4 Gauss values; GNU
  uses the bilinear interpolation `1,x,y,xy` (`miscel.cc`
  `sri_recovery_factor()`, table type 0), as documented in
  `integration_quadrature_and_stress_recovery.md` section 2 (identified
  on probes that may have had no `xi*eta` stress mode). Related: C42
  (physical vs natural coordinates). The kappa probe of F4 T6f WU2
  (`test-2026/kappa_ip_state1`) shows the same signature: reactions equal
  the Professional within 1e-4, nodal sigyy/kap a few per cent off.
- **Next step**: an elastic quad4 probe with a non-zero `xi*eta` stress
  mode (e.g. pure bending of a single quad4) to discriminate the two
  operators; then decide together with C42 (cd_method exit code).
- **Severity**: LOW (nodal output only; the solution is unchanged).
  **Effort**: cheap probe, medium fix.

### C63. `post_node_rhside_fixed` of Pro `slide4` differs from the Professional (F4 T6e, 2026-10-11) — RESOLVED 2026-11-05 (by C110)

- **Resolved by C110 (2026-11-05):** `post_node_rhside_fixed_free()` is re-run after the post-solve pass; Pro `other/slide4` now gives `post_node_rhside_fixed` -0.00577060705357 0.00999500249875 = Professional -5.770607e-3 / 9.995002e-3 (= `post_node_result 10`).

- **Re-measured after C62 (2026-10-17)**: unchanged (GNU 3.44044e-3 /
  7.88353e-3, Professional -5.770607e-3 / 9.995002e-3 = its
  `post_node_result 10`, which the GNU reproduces). `post_node_rhside_fixed_free()`
  runs inside the equilibrium iteration on the pre-solve residual, not
  after the post-solve pass. Candidate: call it once more after the pass
  (`top.cc`). DIVERGENCES C63 row added.

- **Symptom**: Pro `other/slide4`: GNU `post_node_rhside_fixed` 3.44e-3 /
  7.88e-3, Professional -5.770607e-3 / 9.995002e-3 (= its
  `post_node_result 10` sum of the top, which GNU reproduces to 1e-9).
  Also different on HEAD before F4 T6e (3.38e-3 / 7.73e-3), so not caused
  by the C53 pass.
- **Cause**: not established (which fixed dofs / which rhside terms the
  sum includes when slide nodes exist).
- **Severity**: LOW. **Effort**: cheap/medium.

### C62. The post-solve pass skips `parallel_new_dof_before`: reported `node_rhside` without node_force / nodal mass terms (F4 T6e, 2026-10-11) — RESOLVED 2026-10-17 (`a2e68b6`)

- **Resolution (F4 T6f LOW queue Q3)**: the pass also lost the
  `bounda_force` loads (`bounda()` assigns them after
  `iteration_start()`; the pass zeroes NODE_RHSIDE and does not re-run
  `bounda()`). `top.cc` restores a snapshot of the iteration's
  post-`bounda()` NODE_RHSIDE and runs `parallel_new_dof_before` in the
  pass. node_dof identical in the whole corpus (omen A/B, 1019 runs);
  node_rhside moved toward the Professional in 62 of the 64 changed cases
  with a Professional run (e.g. `nodfor1`, `truss6`, `interface1` now =
  Professional); `trubea2`/`trubea3` (GNU planar beam, C37, already rc
  1) now report the applied 0.01 at the loaded node. Tests
  `test-2026/node_rhside_post_pass_{node_force,node_mass,bounda_force}1`.

- **Symptom**: an equilibrium iteration runs `element_loop();
  parallel_contact; parallel_new_dof_before; slide(); solve()`; the
  post-solve material state pass (`top.cc`) runs `element_loop();
  parallel_contact; slide(); parallel_new_dof_diagonal` without
  `parallel_new_dof_before`, so the `node_rhside` written by the pass
  (the step output) lacks the `node_force` and nodal mass terms that
  routine adds. Pre-existing since F4 T5 (C35); seen while closing C53.
- **Fix candidate**: call `parallel_new_dof_before` in the pass too
  (check that it does not double-apply anything to the principal dofs),
  then probes + corpus A/B.
- **Severity**: LOW. **Effort**: cheap/medium.

### C61. Interface Mohr-Coulomb under `inertia_apply`: first-step velocity differs from the Professional (F4 T6d, 2026-10-11) — DEFERRED (class B, F4 T6f LOW queue, owner plan 2026-10-17)

- **Deferred (2026-10-17)**: cause not established (interface + inertia first step, medium);
  no corpus case affected. May share the orientation question of C60.

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

### C60. Interface `element_interface_intpnt_stress` normal sign depends on the node order (F4 T6d, 2026-10-11) — RESOLVED 2026-10-19 (owner decision: adopt the Professional rule)

- **Fix (F4 T6f C60, 2026-10-19)**: `interface_orient_normal()`
  (`interface.cc`) points the normal from side 1 towards side 2 by the
  geometry (side-centroid offset of a thick interface plus the attached
  continuum element of each side, `NODE_ELEMENT`); a free zero-thickness
  interface keeps the default normal. Compression is a NEGATIVE normal
  strain for every node order; the tension cut-off now opens on a
  positive strain. Professional 02-08-2026 = GNU on zero-thickness and
  thick interfaces between blocks (four node orders each: -9.999999 /
  +/-2) and on the free interface orientation (-1 under compression for
  `1 2 3 4` / `4 3 2 1`). Test-2026 `interface_orientation_node_order1`
  (RED on 6bc7457: +10 / +10 / 0 / 0; GREEN -10 everywhere). Five
  GNU-authored test-2014 targets written with the old orientation
  corrected (`iface_mc_gap`, `iface_mc_dil`, `iface_mc_dil_1step`,
  `iface_expansion`, `iface_tangref`; details in the developer manual
  `group_interface.md`). Residual: the free interface pulled apart
  carries tension in the GNU and 0 in the Professional (C87).

- **Re-scoped (F4 T6f LOW queue Q2, 2026-10-17, stopped before any code
  change)**: the premise "output-only sign leak" is false. Probes
  (`scratchpad/lowq/c60/`): two quad4 blocks with a quad4 interface,
  right block pushed left. Elastic: Professional normal -9.98 for all four
  node orders (`2 4 9 10`, `9 10 2 4`, `4 2 10 9`, `10 9 4 2`); GNU +9.98,
  +9.98, -9.98, -9.98. With `group_interface_materi_plasti_mohr_coul_direct
  1 0.2 0 0` + `_tension_direct 1 1`: Professional normal -9.98, shear
  2.023 for every order; GNU orders 1-2 = Professional in magnitude
  (sign flipped), orders 3-4 treat the compressed interface as OPEN
  (normal 0, shear 2.027, block stress 4.05 instead of -5.90). The GNU
  internal convention (compression = positive strain, normal from
  the node numbering rule `nodes[ns1] < nodes[0]`) thus changes the
  physics for half the node orders. Free interface without blocks
  (`1 2 3 4`, `3 4 1 2`, `2 1 4 3`, `4 3 2 1`, side 2 pushed down):
  Professional -1 / 0 / 0 / -1 (pure numbering rule; an opened elastic
  interface carries 0), GNU -1 / -1 / +1 / +1. So the Professional
  orients by geometry when continuum elements are attached and by
  numbering otherwise, and an opened elastic interface carries no
  tension. A fix touches the interface orientation (`interface.cc`
  ~210-245) and the tension behaviour: not a cheap item.

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

### C53. Post-solve material state pass for slide models (`slide4`, F4 T5, 2026-10-10) — RESOLVED 2026-10-11 (F4 T6e)

- **Symptom**: with the ungated pass Pro `slide4` (quad4 on a frictional
  `slide_geometry`, `slide_plasti_friction`, `slide_stiffness`) moved
  away from the Professional (velocity 0.49 relative; before 0.024 sig).
- **Root cause** (read-only investigation, A/B on omen `f4t6e_runs`):
  (1) the pass did not re-apply the slide and contact forces, which are
  added outside the element loop (`slide()`, `parallel_contact`);
  (2) slide law: an absent `slide_plasti_tension` defaulted to
  `sig_t = 0` (a compression-only gap, contrary to manual 6.1044 and the
  code comment), and the friction cap was `max(c + mu*S, 0)` while the
  Professional caps with `c + tan(phi)*|Fn|`; (3) `slide_penalty_law`
  indexed `node_lhside`/`node_rhside` with the node_dof index
  `vel_indx+idim*nder` (the y term went into the eptxx slot of
  test-2014 `node_slide`).
- **Fix**: `slide.cc` (sig_t 1.e300, `|S|` cap, principal-dof index
  `vel_indx/nder+idim`); `top.cc` pass calls `parallel_contact` and
  `slide()` after its element loop; `post_solve_state_pass_gated()` and
  its declarations deleted (DIVERGENCES `POSTPASS-GATE` removed).
- **Verification**: Pro `slide4` = Professional (node 1
  `node_slide_force` -2.11291867e-3 / -1.21989198e-3 vs -2.11291887e-3 /
  -1.21989209e-3; HEAD rc 1); new test-2026 `slide_friction_cohesion1`
  (Pro -2.065482222595e-3; HEAD -0.002 rc 1, GREEN -2.06548203e-3) and
  `contact_block_pass1` (Pro node_dof within 7e-13; HEAD rc 1);
  test-2014 `node_slide` recalibrated to the analytic penalty reaction
  -5e13 + eptxx 0 (owner decision, CORPUS-STALE-TARGETS §7). Out of
  scope, recorded: C62, C63.
- **Severity**: LOW. **Effort**: medium.

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
| B15 | sfnet | `truss2`/`truss3` targets calibrated on the legacy `-updated` truss memory default (F4 T6d) — resolved 2026-10-11 (inputs adapted as the Pro copies) | closed | - |
| C1 | sfnet | `pruebavel` dead/malformed | LOW | cheap |
| C2 | validation-suite | `groundflow_interface` `gvelx` sign | LOW | medium |
| C3 | sfnet | `gforce10` 3D section sign — sign resolved 2026-10-11 (F4 T6a: plot_switch convention, test-2014 input restored, -5000 = Pro); square-face tie kept as deliberate divergence (owner 2026-10-11) | closed | - |
| C5 | — | `contact_target_geometry` index-0 alias — resolved 2026-10-07 (F2 T3a) | MEDIUM | medium |
| C6 | validation-suite | `iface_gen*` generation-coverage gap — resolved 2026-10-11 (F4 T6a: generated elements removed from the inputs) | MEDIUM | cheap |
| C7 | — | Mazars damage vs Pro: Pro equivalent strain weights 4,2,1 (also with plasticity, C7b) — closed 2026-10-12 as a deliberate divergence (class A, F4 T6f); shared non-objective gamma noted, no fix | closed | - |
| C8 | — | default `mesh` `-follow_material` vs Pro `-fixed_in_space` — resolved 2026-10-07 (F3 T2) | MEDIUM | large |
| C9 | — | `-total_linear` + compressibility: Pro time-step dependent (deliberate divergence) | LOW | cheap |
| C10 | Pro `other` | static pore pressure tied to `post_calcul -total_pressure` in the Pro (`excavate1` target) | MEDIUM | cheap |
| C11 | — | `control_materi_dynamic` and implicit truss dynamics vs Pro — resolved 2026-10-11 (F4 T6d) | MEDIUM | medium |
| C12 | — | bicg stop does not abort the GNU on a singular system — closed 2026-10-17 (superseded: re-measured, the Professional 02-08-2026 does not stop either; GNU C77 LU fallback with a warning) | closed | - |
| C13 | — | `control_solver_bicg_stop` at the timestep index: Pro reads index 0 only (deliberate divergence) | LOW | cheap |
| C14 | — | `element_normal`: curved quadratic sides, `node_normal`/`element_boundary`/`global_normal` not implemented (zero normal and plane MC direct return resolved) | LOW | cheap |
| C15 | — | plane direct cut-off phi_flow factor M/G vs Pro 2 (deliberate divergence) | LOW | cheap |
| C18 | — | default memory with `materi_displacement`: Pro behaves as `-total_linear` — resolved 2026-10-07 (F3 T2) | MEDIUM | large |
| C19 | — | explicit `-updated_linear` on `-follow_material`: GNU reference gradient vs Pro current — resolved 2026-10-07 (F3 T2) | LOW | medium |
| C20 | — | `-updated_jaumann` not implemented — resolved 2026-10-08 (F3 T2) | LOW | medium |
| C21 | — | `control_materi_updated_apply -yes`: Pro ignores its manual (deliberate divergence) | LOW | cheap |
| C22 | Pro `other` | `hypo12` marginal sigyy target miss after the default-memory change — deferred (class B, 2026-10-17) | LOW | cheap |
| C25 | — | uninitialized-buffer audit follow-ups — `pef_vel_prev` realloc resolved 2026-10-17 (`efff819`), beam/interface already fixed; NODE_EIGEN print, split, mpc, nonloc deferred (class B) | LOW | cheap |
| C26 | — | `control_print_dof -separate_index`: Pro writes one file per dof, GNU one `dof.N` | LOW | cheap/medium |
| C27 | — | FRD/gmsh layout vs Pro (PSTEP per block, `ALL` component, names, extra blocks/fields, multi-component gmsh fields) | LOW | medium |
| C29 | Pro `other` | `group_condif_conductivity` per-direction values (condif11) not parsed | LOW | cheap/medium |
| C30 | Pro `other` | `mohr_coul_direct4`: Pro does not enforce the interface MC strength (deliberate divergence, GNU = manual 6.632) | LOW | — |
| C31 | Pro `other` | `tet10` nodal quadrature only degree 1 (`tet10_2` gravity) — resolved 2026-10-09 (F4 T4: Keast 5-point + LSQ recovery + HRZ lumping) | MEDIUM | medium |
| C32 | Pro `other` | `axisym9` (stale row, already = Pro), `distri2` (random generator): not GNU defects — closed 2026-10-09 (F4 T4) | LOW | — |
| C33 | — | `_visco` direct cut-offs: `DTIME` from the previous step, no visco on the principal tension cap — resolved 2026-10-09 (F4 T3) | MEDIUM | cheap/medium |
| C34 | — | direct cut shear plastic strain: Pro twice the tensor value (deliberate divergence) | LOW | — |
| C35 | — | one-iteration lag of the whole material state (K02, `spring4` 1 iteration) — resolved 2026-10-10 (F4 T5, post-solve material state pass); inertia/slide/contact widened 2026-10-11 (C52, C53) | MEDIUM | medium/large |
| C36 | validation-suite | `young_power` + `young`: Pro adds both stiffnesses, manual says law = modulus (K06) — decided 2026-10-06: deliberate divergence (GNU = manual), no code change | LOW | cheap |
| C37 | Pro `other` | planar GNU beam vs Pro 3D frame beam (`trubea2/3`, K17) — RESOLVED 2026-11-02 (native 3D beam = Professional) | MEDIUM | large |
| C38 | Pro `other` | `control_support_edge_normal_stiffness_freeze` parse-only (support_edge_normal_3, 1e-9 effect) | LOW | cheap/medium |
| C39 | Pro `other` | K17 residuals: `truss8` birth node resolved 2026-10-06 (death velocity reset open), `update2` -updated robustness, `truss13` truss_distribute_mpc, zero-force truss first iteration | LOW | cheap/large |
| C40 | Pro `other` | element-target contact residuals after F2 T3b: no Pro check "target element has all nodes on the boundary", hex8 target sides never detected, (element_boundary / mesh_boundary resolved F5 T2 V3) | LOW | cheap |
| C41 | Pro `other` | `mohr_coul_direct7`: `-total_linear` + velocity `mpc_node_number` tie — tie resolved 2026-10-09 (F4 T2); residual disy/sigzz difference | LOW | medium |
| C42 | Pro `other` | quad4 nodal stress recovery: natural bilinear (GNU) vs physical linear least squares (Pro) — resolved 2026-10-12 (F4 T6f K21: `axisym2` sigzz 10759.2155 = Pro) | closed | - |
| C43 | Pro `other` | memory hogs of the Bi-CG solver (`force11`, `large1`, `ground19`, `mesh_activate1`, `foundation_consolidation_settlement`): resolved 2026-10-09 (F4 T7 A dense n x n matrix, B lazy band retry, C macro per step); solver speed resolved 2026-10-10 (F4 T9: `-matrix_pardiso` -> SuperLU, CSR operator, recovery table); `foundation...` still Bi-CG bound | MEDIUM | cheap/medium |
| C44 | Pro `other` | one-shot mesh records (`control_mesh_macro`/`_copy`/`_mirror`/`_switch`) re-applied on every time step — resolved 2026-10-09 (F4 T7 C) | MEDIUM | cheap |
| C45 | — | `control_mesh_copy`/`_mirror` after the first time step: `node_dof_start_refined` error — resolved 2026-10-10 (F4 T8 `0ca144c`) | LOW | cheap |
| C46 | Pro `other` | `step_start()` mesh residuals — F4 T8: order, move once, switch pair, `mesh_activate1` resolved; rotate sweep parked, not scheduled (C46c, large); move/switch `node_start_refined` deliberate divergences (A, 2026-10-10) | LOW | large (C46c) |
| C47 | Pro `other` | `force11` 3D section moment `mom1y_sig` 5.80223 vs Pro 5.76460 — resolved 2026-10-11 (F4 T12: quad9 -> hex27 extrusion shares the inner-plane nodes) | LOW | cheap |
| C49 | — | build at `-O1` — resolved 2026-10-11 (F4 T13: 34 scratch-buffer call sites fixed, default `-O2`, results identical) | LOW | medium |
| C50 | Pro other | `interface_bar2_hex8` extruded/converted interface mesh numbering | LOW | medium |
| C51 | validation-suite | hardsoil `gamma_p` read the nodal kappa of the previous evaluation (iteration-count oscillation, `mhardsoil_plast`) — resolved 2026-10-10 (F4 T5); cap/cap1/void reads resolved 2026-10-12 (`70235c6`, C68) | closed | - |
| C52 | Pro `other` | post-solve state pass for inertia models (`dynamic3`) — resolved 2026-10-11 (F4 T6d: interface history fix + inertia removed from the gate) | LOW | medium |
| C53 | Pro `other` | post-solve state pass for slide/contact models (`slide4`) — resolved 2026-10-11 (F4 T6e: slide/contact in the pass, gate removed, slide tension default / `\|Fn\|` cap / penalty index) | LOW | medium |
| C54 | — | run-time record names read without a type check into a 1-value dummy (`post.cc`, `print_el.cc`, `repeat.cc`, `data.cc`, `calcul.cc`) — resolved 2026-10-17 (`8edaf3b`: db_type() check, input error; DIVERGENCES C54-TYPE) | closed | - |
| C55 | validation-suite | `print_unknownsrhside()` indexes `node_rhside` (npuknwn values) with `iuknwn` (heap read past the end, wrong column when nder > 1)  — resolved 2026-10-11 (F4 T6a)| LOW | cheap |
| C56 | Pro `other` | input tokens longer than `MCHAR`-1 (`cin >> str`) overflow the parser buffers (`ground19_water_under_dam`, 141-char URL in a comment)  — resolved 2026-10-11 (F4 T6a)| LOW | cheap |
| C57 | Pro `other` | `post_point` on a quad9 -> hex27 extruded mesh — resolved 2026-10-13 (F4 T6f C57: inner-plane `node_start_refined` z; nodal stress part by C70) | closed | - |
| C58 | Pro `validation` | `control_restart`: interface/spring/contact-spring/tendon/damage element history records not verified after the F4 T6c integration-point reset (`number_timesteps` record absent) — interface/spring verified = Pro 2026-10-17 (`ef9af63`); contact-spring/tendon/damage and `number_timesteps` deferred (class B) | LOW | cheap/medium |
| C59 | sfnet | `group_truss_memory -updated` + `materi_displacement` (and `-total*` without it): input error in the Pro — resolved 2026-10-11 (F4 T6d) | LOW | cheap |
| C63 | — | Pro `slide4` `post_node_rhside_fixed` 3.44e-3 / 7.88e-3 vs Pro -5.77e-3 / 9.995e-3 (pre-existing; unchanged by C62) — RESOLVED 2026-11-05 by C110 (= Professional) | LOW | cheap/medium |
| C62 | — | post-solve pass skips `parallel_new_dof_before`: reported `node_rhside` without node_force / nodal mass terms — resolved 2026-10-17 (`a2e68b6`: also the `bounda_force` loads; node_dof unchanged) | closed | - |
| C61 | — | interface Mohr-Coulomb + `inertia_apply`: first-step velocity differs from the Pro (0.412 vs 0.273), independent of the pass — deferred (class B, 2026-10-17) | LOW | medium |
| C60 | — | interface `element_interface_intpnt_stress` normal sign depends on the node order — resolved 2026-10-19 (geometric orientation, side 1 to side 2, = Professional for every node order) | closed | - |
| C88 | Pro probe | spring element elongation / tangent / `-spring1` indices — resolved 2026-10-19 (F4 T6f C88: projection on the fixed direction = Professional) | closed | - |
| C87 | Pro probe | free elastic interface pulled apart: default gap and shear on opened interfaces — resolved 2026-10-21 (F4 T6f C87: default gap 0, no stress on an opened interface = Professional) | closed | - |
| C89 | Pro probe | Mohr-Coulomb interface dilatancy with the gap: GNU differs from the Professional with gap < 10 and several steps — DEFERRED (class B, owner 2026-11-01) | LOW | medium |
| C90 | Pro probe | `control_print_history` file names differ from the Professional | LOW | small |
| C91 | corpus diff | A `-yes` read into a DOUBLE PRECISION record (`force_edge_water`) is stored as a raw enum number | LOW | small |
| C92 | timing | Threads are created per parallel region (small models 5.5x slower on 16 threads); `processors_maximum` default `-no` until F6 — DEFERRED (class B, owner 2026-11-02) | LOW | medium |
| C93 | Pro probe | `post_element_force` is 0 with `global_element_dof_apply -no` (Professional 3.2967) | LOW | small |
| C94 | Pro probe | Automatic time stepping: the Professional step-size law (x1.25, /4) is not reproduced; adjust_minimum_iterations not consumed — DEFERRED (class B, F5 W1 U4) | LOW | medium |
| C95 | Pro probe | `mrange` not checked against the run-time expansion of value ranges (Professional not measured) | LOW | small |
| C96 | Pro probe | Integration-point state of an element added by `control_input` differs from the Professional (CTRL-INPUT-NEWELEM) | LOW | small |
| C97 | Pro example | `tutorial_4` does not run in the GNU (`-force` with `post_point`, GiD print) and its nonlinear result differs | LOW | medium |
| C98 | Pro probe | `print_vtk_coord` / `control_print_vtk_coord -no`: the Professional drops its `coord` point field, the GNU omits `POINTS` — RESOLVED 2026-11-02 (= Professional) | LOW | small |
| C99 | Pro probe | vtk / tecplot file names (`tn<index>.vtk`, `tecplot.plt`) differ from the Professional `<name>.vtk` / `<name>.plt` — RESOLVED 2026-11-02 (= Professional) | LOW | small |
| C100 | Pro probe | `global_element_dof_from_node_dof` default: GNU -yes, manual / Professional -no — RESOLVED 2026-11-03 (default -no = Professional) | LOW | small |
| C101 | Pro probe | `global_node_dof_empty -yes` keeps the displacement (Professional) while the manual zeroes every dof — RESOLVED 2026-11-03 (kept, = Professional) | LOW | small |
| C102 | Pro probe | node only connected to empty elements keeps its velocity in the GNU (Professional 0) | LOW | small |
| C103 | Pro probe | `dof_element_dof`: the Professional also marks work, f, void -yes | LOW | medium |
| C104 | Pro probe | `dof_limit`: prescribed value outside its limit, clamped in the GNU post-solve element pass | LOW | small |
| C105 | Pro probe | `element_intpnt_dof` temperature slot 0 with condif + mechanics (Professional interpolated) | LOW | small |
| C106 | code reading | `control_reset_value_dof` + spatial value record share one buffer (diagram overwritten after the first point) | LOW | small |
| C107 | F5 T2 V7 probes | older per-model plastic strain dofs filled for every model with the h_lump weight (Professional: 0 under another model) — RESOLVED 2026-11-04 | LOW | small |
| C115 | C114 verifier | bounda wall detection broader than the Professional; no wall reduction for hardsoil / camclay / hypo / tension; no per-index factors | LOW | medium |
| C116 | F5 T3 D2 | two plasticity models in one group (tension + vonmises) differ from the Professional | LOW | medium |
| C114 | F5 T3 C2 | wall friction reduction: cohesion not reduced; direct paths ignore / default 1 the factor (manual and Professional: 2/3 on phi, phi_flow, c) — RESOLVED 2026-11-06 | LOW | small |
| C113 | F5 T3 C1 | materi_plasti_hypo_substepping -no ignored by the Masin routes (Professional: small change) — DEFERRED (class B) | LOW | unknown |
| C112 | F5 T3 B | materi_plasti_maximum_iterations default 20000 in the manual, GNU MAX_ITER 1000 — RESOLVED 2026-11-06 (20000, owner) | LOW | small |
| C111 | F5 T2 V9 | post_point_eps_iso / post_point_node_type -node have no effect in the Professional; GNU follows the manual — DECIDED 2026-11-05 (class A, owner) | LOW | small |
| C110 | F5 T2 V9 | post_node_rhside_fixed/_free: GNU largest absolute value, manual and Professional the sum (ratio differs) — RESOLVED 2026-11-05 (sums, Professional ratio rule) | LOW | small (large blast radius) |
| C109 | F5 T2 V8 | Generic post_calcul item names (`aept`) differ from the Professional post_calcul_label (`a_ept`); GNU follows the manual example — DECIDED 2026-11-04 (class A, owner) | LOW | small |
| C108 | F5 T2 C107 probes | V7 per-model dof labels and total strain kappa measures differ from the Professional (the V7 "adds no dof" claim was wrong) — RESOLVED 2026-11-04 | LOW | medium |
| C64 | — | quad4 per-element nodal recovery of the Professional — resolved 2026-10-12 (F4 T6f K21: physical linear fit, 1-point value, on-axis rule) | closed | - |
| C70 | — | 3x3 quad9 / hex8 / hex27 / tet10 / prism6 recovery in physical coordinates — resolved 2026-10-13 (F4 T6f C70) | closed | - |
| C71 | sfnet | `tube1` (`-lobatto` + axisymmetric) NaN — resolved 2026-10-13 (F4 T6f C71: on-axis Lobatto element integrated with Gauss = Pro) | closed | - |
| C79 | sfnet | Lobatto nodal recovery volume weighted (`tube1` r=1.25 column) — resolved 2026-10-13 (F4 T6f C79: point value / corner fit, arithmetic mean = Pro; `msf_element_rule()` on-axis switch) | closed | - |
| C80 | — | `group_integration_method -reduced` refused by the GNU — resolved 2026-10-16 (F4 T6f C80: Pro rule, linear-fit recovery; default factor 10 = manual, class A) | — | — |
| C80.1 | — | selective / B-bar integration, explicit opt-in GNU extension — resolved 2026-10-16 (F4 T6f C80.1: `-volumetric` of `group_element_selective_reduced_integration`; C80.2 added `-lambda` and removed `-both`) | — | — |
| C72 | — | convected (Eulerian) stress: nodal output resolved 2026-10-13 (`c1a1ae9`, absolute recovery = Pro); neither code convects the IP state; remaining Pro-side residuals: Stokes gap, `cd_method` rotation sweep (the `validation_1` drift C72-V1 is class D, decided 2026-10-18) | LOW | medium |
| C73 | — | optional: convect the integration-point state of history materials on an Eulerian mesh (neither code does; class A improvement, owner decision) | LOW | large |
| C74 | — | `validation_1` prints "Bi-CG iteration broke down" about 600 times — resolved 2026-10-13 (F4 T6f C74: one warning per LU-rescued solve; solver part -> C77) | closed | - |
| C77 | — | Bi-CG absolute breakdown test and convergence floor — resolved 2026-10-15 (F4 T6f C77: relative stop, relative breakdown, guarded SuperLU / band LU retry; BiCGStab choice pending) | closed | - |
| C81 | — | ILU / AMG preconditioning of the iterative solver (F6 profiling, owner 2026-10-15) | LOW | medium/large |
| C75 | — | `prism6` integration points: Professional triangle corners, GNU interior (2/3,1/6,1/6) points — decided 2026-10-13: deliberate divergence (GNU half the bending error) | closed | - |
| C76 | — | `prism6` `group_integration_points -minimal` ignored — resolved 2026-10-16 (`dc4a59a`: Professional 2-point rule, `sum h sig / sum h` recovery) | closed | - |
| C78 | — | `post_point` in curved quad9 / hex27: Professional corner-node (sub-parametric) location, GNU isoparametric — decided 2026-10-13: deliberate divergence (class A) | closed | - |
| C65 | — | Professional `post_point_dof` at a node differs from its own `node_dof` — deferred (class B, folded into DIVERGENCES PPD) | LOW | cheap |
| C66 | — | tria6 integration rule — resolved 2026-10-16 (F4 T6f C66 phase 2, `175faf1`: manual 7 / 4 points, linear-fit recovery = Pro, axisymmetric accepted; HRZ lumping, class A TRIA6-HRZ) | closed | - |
| C84 | — | `post_element_force` 2D/3D on tria/tet/prism elements, 2D record 7 vs 6 values and NaN — resolved 2026-10-17 (F4 T6f C84, `f2ae525`: pol() integration, 6-value 2D record, 3-value 2D result, axisymmetric ring forces; non-continuum elements refused, DIVERGENCES C84-REFUSE) | closed | - |
| C85 | — | tria6 h-weighted lumped terms (groundflow velocity, condif absorption) — resolved 2026-10-17 (F4 T6f C85: HRZ weight, TRIA6-HRZ) | closed | - |
| C86 | Pro `other` probe | `post_element_force` with a truss in the selection — resolved 2026-10-19 (F4 T6f C86: trusses and springs included = Professional; beams/contact springs/interfaces still refused, C84-REFUSE) | closed | - |
| C82 | — | `prism15` `-minimal`: Professional 8 points with an indefinite stiffness, GNU 21 — decided 2026-10-16: deliberate divergence (class A) | closed | - |
| C83 | — | tria3/tet10 `-maximal` manual point counts (4/10) vs both codes (3/5) — decided 2026-10-16: deliberate divergence from the manual (class A); axisymmetric tria3 `-maximal` resolved 2026-10-16 (`72f78bd`: 3 nodal points, 1 on the axis, = Professional) | closed | - |
| C67 | — | nodal recovery of the unmeasured element rules — resolved 2026-10-16 (`d875874` bars/quad16, `dc4a59a` tet4 `-maximal`, prism15, prism6 `-minimal`) | closed | - |
| C68 | — | other integration-point state variables (kapsh, cap1pc, rho, damage, Maxwell msig) read back from the nodes — resolved 2026-10-12 (F4 T6f, `70235c6`) | MEDIUM | cheap/medium |
| C68a | — | Professional Maxwell chain differs from the analytic relaxation even homogeneous (msig -23.81 vs -20.25) — decided 2026-10-18: class D, Professional-suspect (guard `maxwell_relaxation_analytic1`) | closed | - |
| C68b | — | cap1 model: Professional c1his -17 vs GNU pc 108, reactions 6x — deferred (class B, 2026-10-17) | LOW | medium/expensive |
| C68c | — | Mohr-Coulomb hardening-softening on kappa_shear 1.2x off the Professional (homogeneous) — deferred (class B, 2026-10-17) | LOW | medium |
| C68d | — | `materi.cc` old rho/epi/epp read with the Maxwell loop index (two or more chains) — resolved 2026-10-12 (F4 T6f, `2327083`) | LOW | cheap |
| C69 | — | plane-stress membrane: yy/zz split of the plastic strain first order in dt (eptzz -0.00721 at dt 0.1 vs analytic -0.00882) — deferred (class B, 2026-10-17) | LOW | medium |

Fifty-three items remain open (F5 T3 D2 added C116, 2026-11-06). Fifty-two items remained open (C115 added by the C114 verifier, 2026-11-06). Fifty-one items remained open (C114 resolved, 2026-11-06). Fifty-two items remained open (F5 T3 C2 added C114, 2026-11-06). Fifty-one items remained open (F5 T3 C1 added C113, deferred as class B, 2026-11-06). Fifty items remained open (C112 resolved, 2026-11-06). Fifty-one items remained open (F5 T3 B added C112, 2026-11-06). Fifty items remained open (C110 and C63 resolved, 2026-11-05). Fifty-two items remained open (C111 decided as a class A divergence, owner 2026-11-05). Fifty-three items remained open (F5 T2 V9 added C110 and C111, 2026-11-04). Fifty-one items remained open (F5 T2 V8 added C109, decided the same day as a class A divergence, owner 2026-11-04). Fifty-one items remained open (C107 and C108 resolved, 2026-11-04). Fifty-two items remained open (F5 T2 V7 added C107, 2026-11-03). Fifty-one items remained open (C100 and C101 resolved, C106 added, 2026-11-03). Fifty-two items remained open (F5 T2 V5b, 2026-11-03, added C100-C105, all LOW and active; the earlier recount follows). Forty-six items remained open (recount of the summary table, 2026-11-02, after C60, C86, C88 and C87 resolved, C89 added and deferred, C90 and C91 added, C92 added and deferred, C93 added, C94 added and deferred, C95, C96 and C97 added, C98 and C99 added and resolved, C37 resolved): one MEDIUM (B1) and forty-five LOW. Thirty-four are active (the sfnet input items B1-B11, B13, B14; C1, C2, C14, C26, C27, C29, C38, C39, C40, C41, C46, C50, C72, C73, C81, C90, C91, C93, C95, C96 and C97; C14, C39, C41, C46 and C72 count because of their open residuals) and twelve are deferred as class B (C22, C25 rest, C58 rest, C61, C63, C65, C68b, C68c, C69, C89, C92, C94). The twelve class A deliberate divergences of the table (C7, C9, C10, C13, C15, C21, C30, C34, C36, C75, C78, C82) are decided, not open; C68a is decided as class D (Professional-suspect, owner 2026-10-18), not open, and so is the C72 residual C72-V1 (C72 stays open for its Stokes and `cd_method` residuals); fifty-two items are resolved or closed. The earlier running count (45 after C84) also counted decided and closed items; the F4 T6f LOW queue resolved C54 and C62, closed C12, partly resolved C25 and C58, re-scoped C60, deferred nine items and added C86. History of the earlier counts: C23 and C24 resolved 2026-10-05; C25, C26, C27 and B10 added; F1 T3 / K04 2026-10-06 added B11, C28-C32; F1 T6 resolved C28, added and resolved B12; F2 T1 added C36; F2 T2 added C37, C38, C39; F2 T3 resolved C5 and added C40; F3 T2 resolved C8, C18, C19, C20 and added B13, C41; F4 T2 fixed the C41 tie and left a residual; F4 T3 resolved C33; F4 T4 closed C32, added C42 and resolved C31; F4 T7-T8 added C43-C46; F4 T9 resolved the C43 solver speed, decided C46b/C46e and added C47, C48, C49; F4 T10 resolved C48 and added C50; F4 T5 added C52, C53 and resolved C51; F4 T13 resolved C49 and added C54-C56; F4 T12 resolved C47 and added C57; F4 T6a resolved C6, C55, C56 and the C3 sign, the C3 square-face tie becoming a deliberate divergence (owner 2026-10-11), and added B14; F4 T6c fixed the `control_restart` integration-point reset and added C58; F4 T6d resolved C11 and C59, fixed the interface history and resolved C52, and added B15, C60 and C61; B15 resolved the same day; F4 T6e resolved C53 and added C62 and C63; F4 T6f added C64-C68 and C69, and closed the K15 `mmc_tension` vertex sigzz/eppzz and `validation_3` eptzz items as deliberate divergences MC-VERTEX-ZZ / MEMBRANE-EPTZZ, owner 2026-10-12; F4 T6f resolved C68 and added C68a-C68d; C68d resolved the same day; F4 T6f closed C7 (Mazars equivalent strain, with C7b) and the AUDIT K14b `hypo3` sigxx reversal as deliberate divergences C7 / K14b, owner 2026-10-12; the `hypo3` hyhis slot residual stays open under AUDIT K14; F4 T6f K21 WU2 widened C72 with the Stokes-flow residual and closed its `cmat_gate` rewrite, no new item; F4 T6f C72 resolved the convected nodal stress output (C72 stays open for the Professional-side residuals) and added C73 and C74; F4 T6f C70 resolved C70 and added C75; F4 T6f C75 closed C75 as a deliberate divergence, owner 2026-10-13, and added C76; F4 T6f C74 resolved C74 and added C77; F4 T6f C57 resolved C57 and added C78, decided the same day as a deliberate divergence, owner 2026-10-13; F4 T6f C71 resolved C71 and added C79; F4 T6f C79 resolved C79 and added C80; F4 T6f C77 resolved C77 and added C81 (F6); F4 T6f C80 resolved C80 and added C80.1; F4 T6f C80.1 resolved C80.1; F4 T6f C80.2 extended C80.1 with `-lambda`, no new item; F4 T6f C67 resolved the bar and quad16 rules and narrowed C67, no new item; F4 T6f C67 remainder + C76 resolved C67 and C76, re-scoped C66 and added C82 and C83; F4 T6f C82 + C83 decided C82 and C83 part 1 as deliberate divergences (owner 2026-10-16) and resolved C83 part 2; F4 T6f C66 measured the Professional tria6 lumping and added C84; F4 T6f C66 phase 2 resolved C66 (lumping class A TRIA6-HRZ, owner 2026-10-16) and added C85; F4 T6f C85 resolved C85; F4 T6f C84 resolved C84 (the unmeasured Professional handling of non-continuum elements is the DIVERGENCES class B row C84-REFUSE); F4 T6f C60 resolved C60 (owner 2026-10-19: geometric interface orientation) and added C87; F4 T6f C86 resolved C86 for trusses and springs (beams, contact springs and interfaces stay in C84-REFUSE) and added C88; F4 T6f C88 resolved C88; F4 T6f C87 resolved C87 and added C89; C89 deferred as class B at the T6f closure, owner 2026-11-01; F5 W1 U3 added C90 and C91; F5 W1 U2 added and deferred C92 and added C93; F5 W1 U4 added and deferred C94; F5 W1 U6 added C95, C96 and C97; F5 W2 V2 added C98 and C99; C98 and C99 resolved the same day (owner: follow the Professional); F5 T2 V4d resolved C37 (owner: native 3D beam). No open issue is CRITICAL or
HIGH. The validation-suite corpus has no open failure — its 233 real tests
pass (`materi_direct_mc` input corrected, C17 resolved). Every remaining failure is a stale or malformed test, an out-of-scope model
or infrastructure limit, a narrow post-processing divergence, a corpus-coverage
gap, or a performance cost. The source defects found during the cycle were fixed
and validated against the corpora. The remaining post-convergence improvements
(items b–f7) are proposals in `IMPROVEMENT-ANALYSIS.md`; item (a), the post-solve
stress pass hoist, is the only one already implemented.
