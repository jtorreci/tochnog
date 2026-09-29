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
| C1 | sfnet | `pruebavel` dead/malformed | LOW | cheap |
| C2 | validation-suite | `groundflow_interface` `gvelx` sign | LOW | medium |
| C3 | sfnet | `gforce10` 3D section sign | MEDIUM | medium |
| C5 | — | `contact_target_geometry` index-0 alias | MEDIUM | medium |
| C6 | validation-suite | `iface_gen*` generation-coverage gap | MEDIUM | cheap |
| C7 | — | Mazars gamma with shear; damage + early plasticity path | MEDIUM | medium |

Fifteen items remain open: five MEDIUM and ten LOW. No open issue is CRITICAL or
HIGH. The validation-suite corpus has no open failure — its 233 real tests all
pass. Every remaining failure is a stale or malformed test, an out-of-scope model
or infrastructure limit, a narrow post-processing divergence, a corpus-coverage
gap, or a performance cost. The source defects found during the cycle were fixed
and validated against the corpora. The remaining post-convergence improvements
(items b–f7) are proposals in `IMPROVEMENT-ANALYSIS.md`; item (a), the post-solve
stress pass hoist, is the only one already implemented.
