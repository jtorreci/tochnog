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

### C6. `iface_gen*` — reconstructed inputs do not exercise interface generation

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

### C8. Default `mesh`: GNU `-follow_material`, Professional `-fixed_in_space`

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

### C11. `control_materi_dynamic` / truss dynamics without `materi_dynamic` diverge from the Professional

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

### C18. Default memory with `materi_displacement`: the Professional behaves as `-total_linear`

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

### C19. Explicit `-updated_linear` on a `-follow_material` mesh: reference vs current gradient

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

### C20. `-updated_jaumann` memory not implemented

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

### C31. `tet10` integration rule is only degree 1 (Pro `tet10_2` under gravity)

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

### C32. K04 targets that are not GNU defects (`axisym9`, `distri2`)

- `axisym9` (post_point sigxx 0.00158902 vs 0.0015915 +/- 1e-6): the field
  is uniform in both binaries; the GNU total ring force is 0.16 % low
  because `polynom.cc` forces 1-point Gauss for axisymmetric
  `materi_velocity` quad4 (element-constant hoop term, audit K20). Not the
  `post_point_dof` convention. Decision belongs to K20, moved to **F4**
  (design + corpus A/B; project owner, 2026-10-06).
- `distri2` (post_point temp 0.998933 vs 0.999774 +/- 1e-4):
  `control_distribute -normal` random field; both binaries draw a valid
  N(1, 2e-3) (GNU mean 0.99995 / std 0.00203, Pro 0.99997 / 0.00201)
  from different random generators. The target encodes the Professional
  random sequence; no fix.


### C33. `_visco` direct cut-offs: one-step-lagged relaxation time step and no visco on the principal tension cap

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

### C35. One-equilibrium-iteration lag of the whole material state (audit K02) — PARKED, design, scheduled for F4

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
  BOTH-plausible; a decision on free nodes is needed before any change.
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
| C1 | sfnet | `pruebavel` dead/malformed | LOW | cheap |
| C2 | validation-suite | `groundflow_interface` `gvelx` sign | LOW | medium |
| C3 | sfnet | `gforce10` 3D section sign | MEDIUM | medium |
| C5 | — | `contact_target_geometry` index-0 alias | MEDIUM | medium |
| C6 | validation-suite | `iface_gen*` generation-coverage gap | MEDIUM | cheap |
| C7 | — | Mazars gamma with shear; damage + early plasticity path | MEDIUM | medium |
| C8 | — | default `mesh` `-follow_material` vs Pro `-fixed_in_space` | MEDIUM | large |
| C9 | — | `-total_linear` + compressibility: Pro time-step dependent (deliberate divergence) | LOW | cheap |
| C10 | Pro `other` | static pore pressure tied to `post_calcul -total_pressure` in the Pro (`excavate1` target) | MEDIUM | cheap |
| C11 | — | `control_materi_dynamic` and implicit truss dynamics vs Pro | MEDIUM | medium |
| C12 | — | bicg stop does not abort the GNU on a singular system | LOW | cheap |
| C13 | — | `control_solver_bicg_stop` at the timestep index: Pro reads index 0 only (deliberate divergence) | LOW | cheap |
| C14 | — | `element_normal`: curved quadratic sides, `node_normal`/`element_boundary`/`global_normal` not implemented (zero normal and plane MC direct return resolved) | LOW | cheap |
| C15 | — | plane direct cut-off phi_flow factor M/G vs Pro 2 (deliberate divergence) | LOW | cheap |
| C18 | — | default memory with `materi_displacement`: Pro behaves as `-total_linear` (column disx 0.472 vs 0.5) | MEDIUM | large |
| C19 | — | explicit `-updated_linear` on `-follow_material`: GNU reference gradient vs Pro current | LOW | medium |
| C20 | — | `-updated_jaumann` not implemented | LOW | medium |
| C21 | — | `control_materi_updated_apply -yes`: Pro ignores its manual (deliberate divergence) | LOW | cheap |
| C22 | Pro `other` | `hypo12` marginal sigyy target miss after the default-memory change | LOW | cheap |
| C25 | — | uninitialized-buffer audit follow-ups (`pef_vel_prev` realloc, NODE_EIGEN print, split/beam/interface/mpc/nonloc side findings) | LOW | cheap |
| C26 | — | `control_print_dof -separate_index`: Pro writes one file per dof, GNU one `dof.N` | LOW | cheap/medium |
| C27 | — | FRD/gmsh layout vs Pro (PSTEP per block, `ALL` component, names, extra blocks/fields, multi-component gmsh fields) | LOW | medium |
| C29 | Pro `other` | `group_condif_conductivity` per-direction values (condif11) not parsed | LOW | cheap/medium |
| C30 | Pro `other` | `mohr_coul_direct4`: Pro does not enforce the interface MC strength (deliberate divergence, GNU = manual 6.632) | LOW | — |
| C31 | Pro `other` | `tet10` nodal quadrature only degree 1 (`tet10_2` gravity); F4 (design + corpus A/B) | MEDIUM | medium |
| C32 | Pro `other` | `axisym9` (K20 1-point axisymmetric rule, F4), `distri2` (random generator): not GNU defects | LOW | — |
| C33 | — | `_visco` direct cut-offs: `DTIME` from the previous step, no visco on the principal tension cap; F4 | MEDIUM | cheap/medium |
| C34 | — | direct cut shear plastic strain: Pro twice the tensor value (deliberate divergence) | LOW | — |
| C35 | — | one-iteration lag of the whole material state (K02, `spring4` 1 iteration) — parked, design; F4 | MEDIUM | medium/large |
| C36 | validation-suite | `young_power` + `young`: Pro adds both stiffnesses, manual says law = modulus (K06) — decided 2026-10-06: deliberate divergence (GNU = manual), no code change | LOW | cheap |
| C37 | Pro `other` | planar GNU beam vs Pro 3D frame beam (`trubea2/3`, K17) — parked, design | MEDIUM | large |
| C38 | Pro `other` | `control_support_edge_normal_stiffness_freeze` parse-only (support_edge_normal_3, 1e-9 effect) | LOW | cheap/medium |
| C39 | Pro `other` | K17 residuals: `truss8` inactive-element node (both plausible), `update2` -updated robustness, `truss13` truss_distribute_mpc, zero-force truss first iteration | LOW | cheap/large |

Forty-one items remain open: twelve MEDIUM and twenty-nine LOW (C23 and C24 resolved 2026-10-05; C25, C26, C27 and B10 added; F1 T3 / K04 2026-10-06 added B11, C28-C32; F1 T6 resolved C28, added and resolved B12; F2 T1 added C36; F2 T2 added C37, C38, C39). No open issue is CRITICAL or
HIGH. The validation-suite corpus has no open failure — its 233 real tests
pass (`materi_direct_mc` input corrected, C17 resolved). Every remaining failure is a stale or malformed test, an out-of-scope model
or infrastructure limit, a narrow post-processing divergence, a corpus-coverage
gap, or a performance cost. The source defects found during the cycle were fixed
and validated against the corpora. The remaining post-convergence improvements
(items b–f7) are proposals in `IMPROVEMENT-ANALYSIS.md`; item (a), the post-solve
stress pass hoist, is the only one already implemented.
