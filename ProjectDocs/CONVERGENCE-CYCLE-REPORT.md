# GNU → Professional convergence cycle — recapitulation

Date: 2026-09-16
Branch: `dev/masin-kernel-surface` @ `e9c40c9`
Corpus basis: sfnet 197 (`external-downloads/sfnet/extracted/test`, untracked)
and validation-suite 233 (`validation-suite/test-2014`, untracked); internal
suite 16/16.

This document is the wrap-up of the GNU→Professional convergence cycle. It
records the goal, the method, the closed items with their commits and
acceptance evidence, the metrics before and after, the mechanisms discovered,
and the hypotheses that were refuted along the way. It is a recapitulation, not
a specification; every claim is grounded in the cycle's documents
(`SEGUIMIENTO-CONVERGENCIA.md`, the `CORPUS-*` and `DIAG-*` notes, and the
`manual-developer/*` pages) and in the Engram session record.

---

## 1. Goal

The project is a GNU fork of Tochnog (sfnet 2014). The Professional version
(F. Lorenzo / Roddeman) is the behavioural reference. The cycle's goal was to
converge the GNU fork's *behaviour* on the two reproducible test corpora —
the in-repo sfnet 197 and the reconstructed validation-suite 233 — by:

1. closing the corpus tests that the GNU fork fails, and
2. for every remaining failure, establishing an evidence-backed classification
   (source bug, deliberate divergence, stale input/target, out-of-scope) rather
   than leaving it unattributed.

Convergence here is measured by the corpus classification (`PASS` = `rc=0` with
the test's own `target_value` met) plus the internal `build_safe.sh` suite.
It is a behavioural-convergence metric, not a claim of bit-identity.

---

## 2. Method

The cycle ran as a sequence of narrow, ordered work units. The invariant rules:

- **Evidence before change.** Every candidate fix started from a measured
  symptom and a measured root cause. Source changes were made only after the
  cause was reproduced, ideally against both the Professional binary
  (02-08-2026, and 25-10-2023 for older paths) and the reference GNU 2014 build
  (F. Lorenzo `Sources.zip`).
- **Falsify before change.** A hypothesis was treated as refuted unless the
  measurement supported it. Several scheduled work units were *not*
  implemented because the evidence refuted them (see §6).
- **Per-family, not per-test.** Fixes were grouped by root cause (a parser
  family, a solver path, a quadrature rule), then validated by a full corpus
  A/B so that the per-test classification diff could be shown to be exactly the
  intended set with zero regressions.
- **Clean builds only.** All measurements used clean `build_safe.sh` builds
  (`--clean` where the filesystem required it) with `LC_ALL=C` and a 45 s
  per-test corpus timeout; the scratch worktrees were kept out of the main
  tree.
- **Documentation in the same work unit.** Each change updated
  `SEGUIMIENTO-CONVERGENCIA.md` and the relevant
  `manual-user/`+`manual-developer/` pages (the project's mandatory
  convention).

The two corpora are the measurement basis. The frequently cited "363" corpus
lives only in ephemeral `/tmp` and is not reproducible in-repo; "~207" was a
PASS count, not a corpus size. The honest, reproducible denominators are
**197** and **233**.

---

## 3. Metrics

| corpus | start of cycle | end of cycle | net |
|---|---|---|---|
| sfnet 197 | 160 PASS / 37 RUNFAIL / 0 PARSE | **187 PASS / 9 RUNFAIL / 1 PARSE** | +27 PASS |
| validation-suite 233 | 201 PASS / 32 RUNFAIL / 0 PARSE¹ | **233 PASS / 0 RUNFAIL / 0 PARSE** | +32 PASS |
| combined 430 | 361 PASS (84.0%)² | **420 PASS (97.7%)** | +59 PASS |
| internal suite | 16/16 | **16/16** | — |

¹ The validation-suite start figure is the first full classification after the
corpus was reconstructed from the archived `.dbs` dumps (the `.dat` inputs had
been deleted; see §5.1). Before reconstruction the corpus was effectively 0
runnable tests.

² Sum of the per-corpus start figures (160 and 201). An earlier draft of this
table showed 388 (90.2%) for the combined start, which equals the sfnet *end*
(187) plus the validation start (201); the per-corpus starts sum to 361. The
combined end of 420/430 (97.7%) is consistent with both corpora's ends
(187 + 233).

Classification caveat (honest): `PASS` is `rc=0` with the test's own
`target_value` met, which the binary enforces in `exit_tn`
(`miscel.cc`). For the reconstructed tests, `rc=0` does not imply bit-identity
with the pre-incident archive: the archive-era binary predates many later
physics commits, so a faithful reconstruction can still drift numerically. The
final sfnet bucket is `9 RUNFAIL / 1 PARSE` rather than `10 RUNFAIL / 0` only
because `pruebavel` is classified as `PARSE` (unknown-record error); the PASS
count is unchanged. The validation-suite directory holds one generated dump
file (`print_mesh_dof.dat`) that is not a test input; the 233 real tests are all
`PASS` and the dump is excluded from the denominator.

---

## 4. Closed items

Each row: the item, its commit(s), and the acceptance evidence. Commits are on
`dev/masin-kernel-surface` (ancestors of `e9c40c9`).

| # | item | commit(s) | acceptance evidence |
|---|---|---|---|
| 1 | **Masin-basic 2005 rescue** — separate basic kernel (hypo7/8/9), Eq. (15) m-tensor, `ts3` cross-term coefficient, Eq. (21) structure evolution + the Mašín-2007 Eq. (8) pyknotropy modification | `d25906f`, `5a98958`, `c945833` | hypo7 `sigyy=-231.789` (Pro -231.809), hypo9 `-231.789` (Pro -231.811), hypo8 `-230.534` / `s=1.609` (Pro -230.588; target -230.58 ±0.1); all rc=0; hypo1-4 `.dbs` byte-identical (gated) |
| 2 | **dynamic8** — `bounda_sine` layout (end_time), `*_sine` consumers, per-operation quadrature, Professional nodal recovery, scale-relative Bi-CG breakdown test | `6888875`, `814c9ee`, `06e175d`, `abf9c29` | final `disx20 = -9.98525437506e-4` vs Pro `-9.985253775415e-4` (abs 6.0e-11); prescribed sine matches the Pro series to 3.3e-11 |
| 3 | **mpc5** — deleted-element matrix/`lhside` zeroing (energy-consistent released load) + absolute nodal-stress recovery + post-solve stress sync | `7687e0a`, `fa71193`, `0b23957`, `dae06a3` | `post_point_dof 10 -sigyy = 0.5` at k=1,2,3,5 (target 0.5 ±1e-3); full `node_dof` field = Pro to ~5e-17 |
| 4 | **delete1-4 + `control_mesh_delete_geometry_stop` family** — full-matrix delete-factor scaling, frozen deleted-element IP state, stop/`_stop_geometry` keywords | `7687e0a`, `0735028`, `a548415` | delete2 rc=0 node 2 `sigxx=0.8` (Pro 0.8); delete3 rc=0 `sigyy=0.600002285218` (target 0.6 ±1e-4); family field = Pro |
| 5 | **validation_2** (Professional's own suite) — coupled 1D consolidation | `6eb399d` (ancestor) | `post_point_dof 1 -pres = -0.634519741323` vs Pro `-0.6345197413352` (10 digits); target -0.63 ±0.05 |
| 6 | **validation_8** — `-veln` normal-velocity `bounda_dof` with generated mpc ties | `f244ac5` (ancestor) | 50/50 `mpc_node_number` + 8/8 `mpc_node_factor` records byte-identical to the Pro; rc=0 |
| 7 | **Per-operation quadrature + Professional nodal stress recovery** | `06e175d` | dynamic8 6e-11; `beam2d_1` corner `sigxx=59.5766` (Pro 59.5710); `force5` `-9.8071801515` (Pro -9.80718); `post4` moment 0.666667 |
| 8 | **`-total_linear` incremental formulation** (manual Pro 2.2.1) — reset survives the solve | `f045d04` | reset-only probe `sigyy=1` = Pro; `total3` (the only sfnet `-total_linear` test) rc=0 `disx=0.5` |
| 9 | **Parse/arity family** — variable-length `control_timestep_iterations_automatic` and `hypo_intergranularstrain` | `8d70054` | 7 tests closed (family); full corpus A/B |
| 10 | **`node_dof` family** — spurious `node_dof` requests in `top.cc` parse and the `so.cc` eigen path | `ddc9a5b` | 6 tests closed (examp7/refine4/examp15, eigen1/eigen3/ho_heat1); 167→173 PASS, 0 regressions |
| 11 | **groundflow velocity sign** — restore `v = +k grad(h)` (manual Pro 8.2) | `0838457` | `condif5` rc=0 `node_dof 2 temp = 1.0` (target 1.0); full node vector = Pro (`gvelx=+10`) |
| 12 | **`refine_l` BAR2 fall-through** — missing `else if` in the element-shape dispatch | `0d0ad0c` | `refine5` rc=0 `post_point_dof 0 -temp = 1` (target 1 ±5e-2); the 2014 reference aborts at the same assert |
| 13 | **Bi-CG direct-LU retry band width** — `diff*nprinc` → `(1+diff)*nprinc` | `6eb399d` | `examp2` `-0.634519741323` = Pro (was `-0.7206`); `validation_2` closed |
| 14 | **Corpus reconstruction from archived `.dbs`** — 192 `.dat` recovered | `e983be8`, `5c80b60` | 148 recovered first pass, +doc-grounded recoveries; derivation validated 11/11 known pairs, 72/72 IP counts |
| 15 | **Stale sfnet input/target alignment** — 8 `inertia_apply -yes`, `matrix8` sign, `incnav5` target | `b83d907` (docs) | 174→184 PASS, 10 RUNFAIL→PASS, 0 regressions |
| 16 | **Interface stiffness guard + reset scoping/timing** — guard the optional `group_interface_materi_elasti_stiffness`; reset at its own index; reset after the state update | `f142bdf`, `a548415` | `iface_condif`, `groundflow_interface`, `creset_geom` closed; delete3 reconciled |
| 17 | **`group_materi_plasti_camclay` variable length 3..4** — the legacy 4-value form was unreachable | `91676ab` | `mc_pressure_min`/`_off` rc=0 exact doc targets; 4-value form now reachable; corpus classification byte-identical |
| 18 | **delete3 regression reconciliation** — refine the `control_reset_dof` per-index gate (one-shot vs step-scoped) | `a548415` | delete3 rc=0 `0.600002285218`; `creset_geom` preserved; sfnet A/B byte-identical |
| 19 | **Absolute nodal-stress recovery** — the Pro assigns the weighted average, GNU accumulated increments | `0b23957` | delete2 `0.8`; delete family + mpc1-6 field = Pro; all 191 sfnet `.dbs` dumps byte-identical |
| 20 | **Stale validation-suite target recalibration + `gforce7q4` rebuild** | `c592843` | 14 targets recalibrated + `gforce7q4` rebuilt as the faithful 2-quad4 Pro model; 201→216 PASS, 15 RUNFAIL→PASS, 0 regressions |
| 21 | **`inertia_apply -yes` reconstruction of the 15 solver-divergence tests** | `c2f7f6c` | 13/15 converge to their stored targets; validation-suite 216→229 PASS; the 2 residuals carry separate causes (BACKLOG-OPEN-ISSUES.md) |
| 22 | **Reconstruction residuals** — `iface_tangref` target recalibrated to the honest per-IP value (`4.9998`); `iface_gen_method_gen` restored to the archive's input/generated element split | `fbcf144` | 2 RUNFAIL → PASS; validation-suite 229 → 231 PASS, 0 regressions; sfnet byte-identical |
| 23 | **Post-solve stress pass hoisted** out of the equilibrium loop (performance; improvement item (a)) | `3957809` (+ `e0a3997` docs) | 1.12–1.25× on `large3`/`dynamic8` (large3 1.23–1.25×, dynamic8 1.12–1.21×, delete3 ~1.14×; small tests startup-dominated); mpc5/delete2/delete3 targets unchanged; corpus classification diff 0; internal 16/16 |
| 24 | **`-follow_material` reference geometry + initial `element_truss_force` self-stress** (`elmethod`, `numit`) | `c8ee32e` | validation-suite 231 → 233 real tests **233/233**; reconstructed inputs byte-identical to the archived `.dbs`; sfnet unchanged |

---

## 5. Key mechanisms discovered

### 5.1 The `.dbs` reconstruction method

The validation-suite `.dat` inputs were lost to a too-wide `rm -f *.dat`
(2026-08-29) and the directory is gitignored. They were recovered from the
archived `.dbs` database dumps, which the binary writes at the end of a run and
which mix the parsed input records (canonicalised, comment-free, reordered) with
run-produced records.

The method (full detail in `CORPUS-RECONSTRUCTION.md`):

- **Keep the input records, drop the derived ones.** 103 record names are
  dropped; the core is the 45 names that appear in the 41 known `.dat`/`.dbs`
  pairs in *no* `.dat` (state/results: `node_dof*`, `element_dof*`,
  `node_rhside*`, `element_matrix*`, `dof_*`, `time_*`, `global_solver_*`,
  `post_*`, ...), extended with the run-produced records seen in the wider dump
  set.
- **Initial-state recovery.** `node_dof_start_refined` holds the input
  `node_dof`; it is re-emitted as `node_dof`.
- **Derive, do not invent, the lost scalars.** `number_of_integration_points`
  is never written by the dump, but it is recoverable from the archived
  `element_dof` and `dof_label` record lengths:

  ```
  nuknwn = len(dof_label) + 1
  number_of_integration_points = len(element_dof) / nuknwn
  ```

  The derivation was cross-checked on the known pairs (11/11 exact) and with a
  fresh K-probe (72/72 agreement).

A critical corollary: **`inertia_apply` is never echoed in a `.dbs`** (0 of 430
dumps). Its absence from an archived dump is therefore not evidence that the
input lacked it. This directly enabled the recovery of the 15 solver-divergence
tests (§5.5).

### 5.2 Per-operation quadrature

The Professional does not use a single global quadrature. The quadrature is a
property of the element integral, and the nodal stress recovery is a separate
operator:

| operation | Professional | GNU (aligned) |
|---|---|---|
| element stiffness, internal force, body load (quad4/hex8) | 2×2 / 2×2×2 Gauss | Gauss |
| element stiffness (quad9/hex27) | 3×3 / 3×3×3 Gauss | Gauss |
| mass matrix | row-sum lumped (quadrature-invariant) | unchanged |
| edge / surface load, interfaces | boundary integration; interfaces keep nodal (Lobatto) | unchanged |
| stress recovery (IP → node) | per-element least-squares fit + arithmetic mean | Professional operator |

The historical GNU default was Lobatto for the element integral plus an
`h`-weighted lumped recovery, which conflated the two. Manual 6.623's
"inertia ⇒ lobatto" clause was measured to describe the **mass** (which is
row-sum lumped and quadrature-invariant), not the stiffness; the 02-08-2026
binary contradicts it for the stiffness. The alignment is in `06e175d`.

### 5.3 The Professional's nodal stress recovery

The `.dbs` dumps expose the operator: `element_intpnt_dof` → `element_dof` →
`node_dof`, with `node_lhside` as the assembled weight. Measured operators:

- quad4: bilinear interpolation `1, x, y, xy` (residual vs `element_dof` 1.4e-7);
- quad9: uniform-weight linear least squares `1, x, y` (1.2e-12);
- hex8: uniform-weight linear least squares `1, x, y, z` (1.0e-8);
- hex27: assumed linear LS (not independently verified).

The operator is fixed in the element natural coordinates and uses **uniform**
weights (the Gauss weights give a different operator). Assembly is the
**arithmetic mean** of the adjacent elements' `element_dof`.

The 2026-09-13 alignment fixed the weights but the GNU still *accumulated*
increments (`node_dof_new += node_rhside/node_lhside`). At a node shared by two
elements whose increments cancel, the reported stress stayed at its step-start
value (delete2 reported 0.6 where the Pro reports 0.8). The Professional
**assigns** the absolute weighted average; the fix (`0b23957`) makes the stress
dofs assign the ratio, with a gate that keeps the historic incremental update
for convected stress dofs and for models with an independent strain unknown.

### 5.4 mpc5 / delete energy consistency

Two distinct mechanisms closed the delete/mpc family:

- **Delete-factor scaling (`7687e0a` → `fa71193`).** A partially deleted element
  stores `element_delete_factor` + a frozen `element_rhside_delete` (the
  released-stress internal force). The Professional stores **no**
  `element_matrix_delete` record. The first fix scaled the full matrix by the
  factor; the final fix recognises that the frozen released load is a
  *constant*, so its tangent is zero: with `ELEMENT_DELETE_FACTOR` active the
  element matrix and `element_lhside` are set to **zero** (not scaled). This
  restores the one-pass Newton cancellation `A = dt·K = dt·BᵀCB` of the
  staggered scheme. The historic `factor·K` left a residual stiffness that
  responded to strain while the internal force was frozen, so the released load
  relaxed geometrically (mpc5 contracted at rate exactly 0.5, reaching the
  target only at ~12 iterations instead of 2).
- **Post-solve stress sync (`dae06a3`).** The `materi_stress` dofs are
  non-principal: they are updated *after* the solve from the *pre-solve*
  constitutive stress, so the reported nodal stress trailed the solved velocity
  by one iteration (mpc5 reported 1.0 at k=1; the Pro reports 0.5). The fix adds
  a stress-only `element_loop()` pass after `solve()` and before
  `parallel_new_dof_diagonal()`, re-evaluating the constitutive stress from the
  solved velocity. This is *not* algebraically the same as iteration k+1: it
  adds a constitutive evaluation, not a second momentum solve.

### 5.5 The `inertia_apply` default gotcha

Commit `0560dcf` aligned the GNU `options_inertia` default with the
Professional: absent `inertia_apply` now means `-no` (manual Pro 6.785), where
the historic GNU default was `-yes`. This deliberate alignment fixed 13 static
corpus tests but silently changed the behaviour of legacy models that relied on
the old default:

- 8 sfnet inputs + `incnav3` + `examp2` needed an explicit `inertia_apply -yes`
  (the Pro's own copies carry it);
- 15 validation-suite models needed it too. Those models have a material
  density but insufficient Dirichlet constraints, so the legacy `-yes` supplied
  a `rho/dtime` lumped mass diagonal that regularised an otherwise **singular**
  static system. Without it the matrix has genuine rigid-body modes, the Bi-CG
  breaks down and the direct-LU retry diverges.

The onset was isolated by a clean-build bisect: `0560dcf` (09-05 07:06) fails
all 15, while its predecessor and the quadrature commit `06e175d` do not. The
`-no` default is a deliberate Professional alignment and was **not** reverted;
the fix is the explicit record in the reconstructed inputs.

### 5.6 The honest solver as an explanation for "divergences"

`so_bicg.cc` previously reported success on three non-convergence exits. The
current honest criteria report the true residual, so several tests that "passed"
on the archive-era binary now correctly report divergence. This is a deliberate
improvement, and it is why some archived targets were met only by a dishonest
success. It also exposed a latent bug: the direct-LU retry band width used
`diff*nprinc` instead of `(1+diff)*nprinc`, so the first-to-last dof coupling of
adjacent nodes wrote outside `mat`/`ldmat` and silently corrupted the retry
(`6eb399d`).

### 5.7 The `-follow_material` reference geometry

Under `-follow_material` the mesh follows the material, so the `node` record in
a `.dbs` dump is the **deformed** geometry, while the input reference is
preserved in `node_start_refined` (`top.cc` copies the input `NODE` there at the
start of the run). The reconstruction had used the deformed `node` as input,
which is why `elmethod` and `numit` would not reproduce their archived targets.
The second lost record was the initial `element_truss_force` self-stress, which
the dump stores only in its relaxed final state (≈0). With the reference
geometry restored and the initial force derived from the documented relaxation
law (`F_init = k·(L_init − L_final)`, with `k = E·A/L_init`), the reconstructed
inputs reproduce the archived `.dbs` byte-for-byte. This is a
corpus-reconstruction mechanism, not a source change.

The decisive evidence was a surviving original `.dat`/`.dbs` pair that uses
`-follow_material`, `dpline_method`: the input `node 3 0 0.98` becomes
`node 3 0 0.96` (deformed) in the dump, while `node_start_refined 3 0 0.98`
preserves the input. The reconstructed `elmethod`/`numit` had used the deformed
`node` (`0.877` / `1.631`) as input; the correct input is the reference
(`0` / `1` / `2`). The initial truss forces are fully determined by the archived
geometry and the relaxation law: elmethod 0.123 / 0.246, numit 0.123.

---

## 6. What made it work — and the refuted hypotheses

The cycle's discipline was to treat every scheduled fix as a hypothesis and to
record the refutations as first-class results. The refuted or corrected
hypotheses were load-bearing:

- **Masin-basic R2 ("reference geometry refuted")** was a *false negative*: the
  parked kernel had a transcription bug (Eq. 15 used `T_hat:T_hat*` instead of
  `T_hat:T_hat`), so the earlier "R2 does not improve the residual" verdict was
  measured with the wrong kernel. With the kernel fixed, R2 is part of the
  hypo7/8/9 closure.
- **"updated_linear log vs linear"** as the cause of the Masin stress residual
  was refuted by the R2 measurement: the strain measure matches the Pro
  (`eptyy = -0.300000000`) but the stress residual does not improve
  systematically.
- **Monolithic saddle-point mixed solve (mpc5 candidate (a))** was refuted as
  structurally inconsistent for shared nodes in the current tie scheme
  (mpc1's shared σ2 cancels to 1.667, not 2.0). Candidate (b) — post-solve
  re-evaluation — was implemented instead.
- **"Candidate (b) is algebraically identical to iteration k+1"** was corrected:
  the *value* sequence coincides, but the *mechanism* (a constitutive
  evaluation, not a solve) is distinct and does sync the k=1 report.
- **"examp6 is pre-existing"** was corrected to a regression from the
  `ELEMENT_DOF` stress-slot fix `8f16dc0`, compounded by `0560dcf`; still
  out-of-scope because the target is tied to a residue-driven adaptive-remesh
  node numbering.
- **"examp24 cost is mesh generation"** was corrected to O(N²) mesh
  **mapping** (`map.cc` `point_el`/`matrix_ab`, 29M `point_el` calls per cycle).
- **validation_2's `-0.688985`** was refuted as a corrupt-retry artifact; the
  honest value is `-0.634519741323` = Pro.
- **validation_8's "not deducible without more data"** was refuted: `-veln` was
  implemented and reproduces the Pro's 50/50 + 8/8 MPC records.
- **The 15 validation-suite solver divergences** were first attributed to the
  honest Bi-CG criteria and the quadrature commit; both were refuted by the
  bisect, and the true cause (`inertia_apply`) was found.
- **The corpus basis itself** was corrected: "363" is ephemeral `/tmp`, "~207"
  was a PASS count; the reproducible denominators are 197 and 233.

The pattern is consistent: the fixes that survived were the ones whose root
cause had been measured against a reference binary, and the refutations
prevented speculative edits that would have regressed the corpora.

---

## 7. Residuals

The remaining non-PASS tests are enumerated and triaged in
`BACKLOG-OPEN-ISSUES.md`; the follow-on opportunities are in
`IMPROVEMENT-ANALYSIS.md`. In summary:

- sfnet 197: 9 RUNFAIL + 1 PARSE (`pruebavel`, dead/malformed). The 9 RUNFAIL
  are all out-of-scope: 4 stale targets (`examp17`, `ho_mech2`, `ho_mech4`,
  `taylor3`), 3 solver/timeout (`examp14`, `examp24`, `ho_mech3`), 1
  contact/refinement (`examp6`) and 1 infrastructure (`petsc1`).
- validation-suite 233: 0 RUNFAIL — all 233 real tests pass. The directory's
  remaining non-PASS file is the generated dump `print_mesh_dof.dat`, not a test
  input. The four residuals open at the cycle snapshot (`iface_tangref`,
  `iface_gen_method_gen`, `elmethod`, `numit`) were closed by `fbcf144` and
  `c8ee32e`; see the resolved section of `BACKLOG-OPEN-ISSUES.md`.

## 8. Sources

- `ProjectDocs/SEGUIMIENTO-CONVERGENCIA.md` — the verification ledger.
- `ProjectDocs/CORPUS-RECONSTRUCTION.md` — the `.dbs` reconstruction method and
  the classification passes.
- `ProjectDocs/CORPUS-STALE-ALIGNMENT.md`, `CORPUS-STALE-TARGETS.md` — the stale
  input/target alignment and the recalibrations.
- `ProjectDocs/CORPUS-SOLVER-TIMEOUT.md`, `CORPUS-PARSE-INIT-EXAMP2.md`,
  `DIAG-PHYSICS-MATERIAL-STATE.md`, `DIAG-CONTACT-REFINEMENT.md` — the triage of
  the sfnet residuals.
- `ProjectDocs/VALIDATION-PROFESSIONAL-2-8.md` — validation_2/8.
- `ProjectDocs/manual-developer/` — the per-mechanism developer pages
  (`integration_quadrature_and_stress_recovery.md`, `nodal_stress_recovery.md`,
  `nodal_stress_post_solve_sync.md`, `mpc5_staggered_stress_lag.md`,
  `delete_factor_element_matrix_scaling.md`,
  `control_mesh_delete_geometry_stop.md`, `solver_bicg_band_width.md`,
  `groundflow_velocity_sign.md`, `refine_locally_bar2.md`,
  `control_reset_dof_scoping.md`, `node_dof_family.md`).
