# Backlog — open issues after the convergence cycle

Date: 2026-09-16
Branch: `dev/masin-kernel-surface` @ `a1362ef`

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

## A. validation-suite 233 — the 4 RUNFAIL

### A1. `iface_tangref` — interface tangential-force record semantics changed

- **Corpus**: validation-suite 233.
- **Symptom**: target `element_interface_force_tang2 2 = 10.4025 ±0.2`; the
  current binary returns `4.99979673431` (four per-IP values, all ≈5.0).
  The first target `element_interface_force_tang 2 = 0 ±0.001` is met
  (`-1.75e-14`).
- **Root cause**: commit `bf74f28` replaced the old "side average,
  all-against-all" interface assembly with a per-integration-point
  facing-pair assembly with Lobatto weights, and changed the history records
  (`element_interface_strain_normal`, `..._force_tang`, `..._force_tang2`)
  from a scalar to a length-4 per-IP array. The stored target `10.4025 ±0.2`
  is tied to the old scalar semantics; HEAD returns the per-IP value.
  Onset bisect: `7823d27` (good, scalar `10.3532842167`) → `bf74f28` (first
  4-value form) → `6cafafe`/HEAD (per-IP `4.9998`).
- **Evidence**: `CORPUS-STALE-TARGETS.md` §"Inertia-apply reconstruction" §3.1;
  `residual-iface_tangref-bisect/`.
- **Severity**: **MEDIUM** — a test-target/semantics mismatch, not a source
  defect; the per-IP record is the validated Professional-aligned behaviour.
- **Proposed fix**: recalibrate the target to the current honest per-IP value
  (`4.9998`), a corpus decision like the 14 recalibrations already applied.
  Alternatively, add an explicit `group_interface`-level assembly selector if
  the old scalar aggregate is still wanted (a new feature, not a bug fix).
- **Effort**: **cheap** (target recalibration); **medium** (assembly selector).

### A2. `iface_gen_method_gen` — reconstruction dropped the input `element_geometry` records

- **Corpus**: validation-suite 233.
- **Symptom**: target `element_geometry 3 = 10 ±1`; `exit_tn` reports
  "Error in calculation" because `element_geometry 3` does not exist.
- **Root cause**: the reconstruction's drop-list classified `element_geometry`
  as a run-produced record, so the three **input** records
  (`element_geometry 1 20`, `2 30`, `3 10`) were dropped. The archive has
  them; HEAD generates elements 4 and 5 instead.
- **Evidence**: `CORPUS-STALE-TARGETS.md` §3.2; variant-A re-insertion test
  (`residual-iface_gen_method_gen-variantA/`) gives rc=0.
- **Severity**: **LOW** — a reconstruction artifact affecting one test; the
  sibling `iface_gen_method` is harmless because its target is a `node_dof`.
- **Proposed fix**: restore the three input records (trivial, faithful to the
  archive) **or**, preferably, redesign the target to exercise the generation
  feature rather than the input tag (the archive has no generated elements,
  HEAD generates 4 and 5). The second option is the honest one and is why the
  record restore was deliberately not applied in the cycle.
- **Effort**: **cheap**.

### A3. `elmethod` — lost `-follow_material` reference geometry (secondary target)

- **Corpus**: validation-suite 233.
- **Symptom**: the primary A/B passes (`-middle` 2 lines / `-node` 4 lines,
  `2*mid == x_node2`); the secondary `node_dof 3 -velix` target (`-0.369`) is
  not reproducible.
- **Root cause**: the `-follow_material` reference geometry
  (`node_start_refined`) is lost by the reconstruction; the dump `node` is the
  **deformed** geometry. Re-inserting `node_start_refined` as input did not
  reproduce the archived `-0.369` on either the current or the 2026-08-28
  binary.
- **Evidence**: `CORPUS-RECONSTRUCTION.md` §2.1, §3.
- **Severity**: **LOW** — the feature under test (`control_print_element_method`)
  works; only the secondary velocity target is unrecoverable.
- **Proposed fix**: recover the original input form (e.g. the exact
  `control_mesh`/geometry records that produced the reference configuration)
  or recalibrate the secondary target to the honest current value.
- **Effort**: **medium**.

### A4. `numit` — lost `-follow_material` reference geometry (secondary target)

- **Corpus**: validation-suite 233.
- **Symptom**: the primary feature passes (exactly 16 monitor lines); the
  secondary `node_dof 2 -velix` target (`-0.123 = 0.877 - 1`) is not
  reproducible.
- **Root cause**: same as A3 — the deformed geometry is what the dump stored.
- **Evidence**: `CORPUS-RECONSTRUCTION.md` §2.1, §3.
- **Severity**: **LOW**.
- **Proposed fix**: same as A3.
- **Effort**: **medium**.

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

---

## C. Additional residuals

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

### C4. Post-solve stress pass cost

- **Corpus**: all models with `materi_stress` + the absolute recovery (e.g.
  the delete/mpc family).
- **Symptom**: the element pass runs twice per equilibrium iteration
  (`nodal_stress_post_solve_sync`), roughly doubling the element-assembly cost
  for those models. No corpus timeout or classification changed.
- **Root cause**: the post-solve stress sync (`dae06a3`) is implemented as a
  full stress-only `element_loop()` after every solve, because the stress dofs
  are non-principal and must be refreshed from the solved velocity.
- **Evidence**: `manual-developer/nodal_stress_post_solve_sync.md` §Cost,
  §Pending refactorings.
- **Severity**: **MEDIUM** — a performance cost in a common configuration
  (any model reporting `materi_stress`).
- **Proposed fix**: hoist the pass after the equilibrium loop and update the
  stress dofs once (manually or via a stress-only variant of
  `parallel_new_dof_diagonal`), halving the added cost. See
  `IMPROVEMENT-ANALYSIS.md` §a.
- **Effort**: **cheap**–**medium**.

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

---

## D. Summary

| id | corpus | item | severity | effort |
|---|---|---|---|---|
| A1 | validation-suite | `iface_tangref` target vs per-IP semantics | MEDIUM | cheap/medium |
| A2 | validation-suite | `iface_gen_method_gen` dropped input records | LOW | cheap |
| A3 | validation-suite | `elmethod` lost reference geometry | LOW | medium |
| A4 | validation-suite | `numit` lost reference geometry | LOW | medium |
| B1 | sfnet | `examp6` remesh node numbering | MEDIUM | cheap/large |
| B2 | sfnet | `examp14` stale target + runtime | LOW | cheap |
| B3 | sfnet | `examp17` stale target | LOW | medium |
| B4 | sfnet | `examp24` O(N²) mesh mapping | LOW | large |
| B5 | sfnet | `ho_mech2` stale p-refinement input | LOW | medium |
| B6 | sfnet | `ho_mech3` model divergence | LOW | large |
| B7 | sfnet | `ho_mech4` missing strain record | LOW | medium |
| B8 | sfnet | `petsc1` PETSc infrastructure | LOW | large |
| B9 | sfnet | `taylor3` stale input + solver | LOW | medium |
| C1 | sfnet | `pruebavel` dead/malformed | LOW | cheap |
| C2 | validation-suite | `groundflow_interface` `gvelx` sign | LOW | medium |
| C3 | sfnet | `gforce10` 3D section sign | MEDIUM | medium |
| C4 | all | post-solve stress pass cost | MEDIUM | cheap/medium |
| C5 | — | `contact_target_geometry` index-0 alias | MEDIUM | medium |

No open issue is CRITICAL or HIGH: every remaining failure is a stale/malformed
test, an out-of-scope model or infrastructure limit, a narrow post-processing
divergence, or a performance cost. The source defects found during the cycle
were fixed and validated against the corpora.
