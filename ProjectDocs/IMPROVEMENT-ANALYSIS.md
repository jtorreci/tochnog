# Improvement analysis — post-convergence

Date: 2026-09-16
Branch: `dev/masin-kernel-surface` @ `a1362ef`

Candidate improvements identified during the GNU→Professional convergence cycle.
Each entry states the problem, the proposed change, the expected benefit, the
risk, and the estimated effort. These are proposals, not commitments; none is
required for the cycle's results.

---

## a. Post-solve stress pass cost — **DONE 2026-09-16** (branch `dev/stress-pass-hoist`)

**Problem.** The post-solve nodal-stress sync (`dae06a3`) runs a stress-only
`element_loop()` after **every** equilibrium iteration for models with
`materi_stress` and the absolute recovery. That roughly doubles the
element-assembly cost for those models (`nodal_stress_post_solve_sync.md`
§Cost). The corpus showed no timeout or classification change, but the cost is
real for production models.

**Proposed change.** Hoist the pass out of the iteration loop: after the
equilibrium loop converges, run the stress-only pass **once** and update the
stress dofs (either directly, or via a stress-only variant of
`parallel_new_dof_diagonal`). The in-loop pass exists only so that any
per-iteration print stays in sync; if per-iteration printing is not required,
the once-per-step form is equivalent for the final report.

**Expected benefit.** Halves the added element-assembly cost for the affected
models, restoring them close to the pre-sync cost.

**Risk.** Low-to-moderate. The pass ordering interacts with `parallel_new_dof_diagonal`
and with the structural-element bypass (`truss`/`beam`/`spring`/`contactspring`/
interface) that requires the non-stress save/restore. Hoisting must preserve
the "only the last report is observed" semantics; a model that prints stress
inside the iteration loop would change output. Mitigation: A/B the full sfnet
197 + validation-suite 233 corpora and the delete/mpc family field.

**Effort.** **cheap–medium** (localised to `top.cc`/`dof.cc`, but needs a careful
A/B).

**Implemented.** Only `top.cc` changed: the save/zero/`element_loop()`/restore
block moved from inside the equilibrium loop to a single block after it, with
the stress dofs assigned directly (stress-only variant of
`parallel_new_dof_diagonal`, including the `materi_stress_pressure_history`
and `derivatives` stress-rate slots). Measured speedups: `large3` (784
elements, 220 iterations) **1.23-1.25x**, `dynamic8` (16 elements, 4000
iterations) **1.12-1.21x**; the small delete/mpc tests are startup-dominated.
Targets unchanged: mpc5 `sigyy=0.5` at `k=1,2,3,5`, delete2 node 2 `sigxx=0.8`,
delete3 `sigyy=0.6`; sfnet 197 and validation-suite 233 classification diff
= 0; internal suite 16/16. Two non-target derived-dof changes are documented
(`contact_heatgen` `xsigyy`, `mpc8` `sigxx`). Detail:
`manual-developer/nodal_stress_post_solve_sync.md` §Cost / §Verification.

---

## b. Solver robustness

**Problem.** The solver is honest but not robust on the genuinely hard systems.
Two distinct effects:

1. **The honest Bi-CG criteria.** `so_bicg.cc` now reports the true residual
   (previously three non-convergence exits were reported as success). This is
   correct, but it means some ill-conditioned systems break down and fall to
   the direct-LU retry, which is the path where the band-width bug lived
   (`6eb399d`). The retry is a correctness safety net, not a scaling solution.
2. **The genuinely hard systems.** The consolidation coupling (`examp2`,
   `validation_2`), the singular static systems that relied on the legacy
   inertia default, and the eulerian convection/plasticity family
   (`validation_1`/`validation_8`'s stress field) are the observed hard cases.

**Proposed change.** Evaluate two directions, in order of cost:

- **Preconditioning.** The Bi-CG has no preconditioner. A simple diagonal
  (Jacobi) or block-Jacobi preconditioner would reduce the iteration count on
  the ill-conditioned systems without changing the equations. The band-width
  fix already removed the silent corruption, so the remaining failures are
  genuine conditioning.
- **MINRES / a mixed solve for the staggered family.** `mpc5_staggered_stress_lag.md`
  candidate (a) proposed a monolithic saddle-point mixed solve
  (`[0 Bᵀ; B −C⁻¹]`) with MINRES or SuperLU with pivoting, making the stress
  dofs principal and removing the one-iteration lag structurally. It was
  **refuted as structurally inconsistent for shared nodes in the current tie
  scheme** (mpc1's shared σ2 cancels to 1.667, not 2.0). It remains the
  principled long-term fix *if* the tie scheme is made energy-consistent for
  shared nodes first; the current energy-consistent tie work (`b3128c0`) is the
  precondition for revisiting it.

**Expected benefit.** Better convergence on the ill-conditioned coupled systems
and, for the mixed-solve direction, a structurally correct stress dof (no
post-solve workaround) and closure of the whole staggered-family residual.

**Risk.** Medium-to-high. A preconditioner changes the iteration path and must
be A/B-validated (the honest criteria make any change visible). The mixed solve
is a corpus-wide refactor of `elem.cc`/`so.cc` with a new solver path and new
validation for the whole suite; the shared-node inconsistency must be solved
first.

**Effort.** **medium** (preconditioner) / **large** (mixed solve).

---

## c. `build_safe.sh` full-mode file-check layer

**Problem.** After the corpus reconstruction, the canary `msf_beam2d.dat` exists
in `validation-suite/test-2014/`, so `build_safe.sh` now runs in **full mode**
and executes the Sprint 11 file-check layer. That layer asserts expected outputs
of files generated during the pre-incident sprints; many of those expected
values are stale relative to the current binary (the post-incident backlog). The
script therefore exits non-zero, and the FALLO list is documented to be
identical to the baseline — i.e. it is not a regression, but it makes the script
a dirty gate that cannot be used as a clean pass/fail signal.

**Proposed change.** Update the stale expected outputs in the file-check layer
to the current honest values (with the same A/B discipline used for the corpus
target recalibrations), or move the file-check expectations into a maintained
golden-file set generated from a known-good binary and reviewed like any other
target. Keep the `rc`-based internal suite (16/16) as the primary gate.

**Expected benefit.** `build_safe.sh` becomes a clean gate again: a non-zero
exit means a real regression, not a stale expectation. This also removes the
need for the "reduced mode" workaround.

**Risk.** Moderate. Blindly recalibrating expectations would mask a real
regression. Mitigation: only recalibrate expectations whose FALLO list is
identical to the baseline; keep the full sfnet/validation-suite corpora as the
independent check; review each expectation's diff before updating.

**Effort.** **medium** (the file-check layer is large, ~1996 lines; the
recalibration is mechanical but must be reviewed).

---

## d. Corpus tooling

**Problem.** The measurement basis is not reproducible from the repository:

- The "363" corpus lives only in ephemeral `/tmp` and is not reproducible
  in-repo; the honest denominators are the sfnet 197
  (`external-downloads/sfnet/extracted/test`) and validation-suite 233
  (`validation-suite/test-2014`), both **untracked/gitignored**.
- The runner is a `/tmp` shell script (`corpus-run.sh` / `run-corpus.sh` /
  `final-corpus-run.sh` variants) with ad-hoc classification output.
- Target recalibration is manual: the current value is measured, then the
  `target_value` line is edited by hand and recorded verbatim in the docs
  (because git does not track the corpus).

**Proposed change.** Three pieces:

1. **Make the corpora reproducible in-repo.** A `corpus/` manifest with
   checksums plus a fetch/extract script (the third-party test data has
   licensing implications, so a fetch script + pinned checksums is safer than
   committing the data). The reconstruction-produced tests are already
   documented and can be regenerated.
2. **A stable runner.** Promote the `/tmp` runner to `scripts/run-corpus.sh`
   with a documented invocation and a machine-readable classification output
   (`name\tPASS|RUNFAIL|PARSE\tvalue`), so results are comparable across runs
   and can feed CI.
3. **Target-recalibration tooling.** A helper that reads a corpus run, extracts
   the observed value for each `target_value`, and proposes the recalibration
   edit with the current tolerance — leaving the human decision (and the
   documented rationale) explicit.

**Expected benefit.** Reproducibility and comparability; the convergence metric
becomes a first-class CI signal instead of a `/tmp` artifact; recalibrations
stop being error-prone manual edits.

**Risk.** Low-to-moderate. The main risk is licensing/redistribution of the
third-party test data; the fetch-script approach avoids it. The runner must
keep the 45 s timeout and the `LC_ALL=C`/`LD_LIBRARY_PATH` discipline to stay
comparable to the recorded results.

**Effort.** **medium**.

---

## e. Reusability of the reconstruction method

**Problem.** The `.dbs` reconstruction method (`CORPUS-RECONSTRUCTION.md`) is
well documented but implemented in ad-hoc `/tmp` scripts (`work/reconstruct.py`,
`work/recon2.py`). If a future incident (or a new dump set) needs the same
recovery, the tooling must be rebuilt from the prose. The method has
version-sensitive parts: the 103-name drop-list and the derived
`number_of_integration_points` formula.

**Proposed change.** Promote the method to `scripts/reconstruct_dbs.py`:
- the drop-list as a data file (name → reason), versioned with the code;
- the `number_of_integration_points` derivation
  (`nuknwn = len(dof_label)+1`, `n_ip = len(element_dof)/nuknwn`) with the
  exact-division check;
- a **validation mode** that runs the derivation against the known
  `.dat`/`.dbs` pairs (the 11/11 exact check) so a regression in the method is
  caught immediately;
- a report of RECOVERED / PARTIAL / UNRECOVERABLE, matching the documented
  classification.

**Expected benefit.** A re-runnable, self-checking recovery tool; the method's
validity is verified on every use rather than asserted in prose. It also makes
the "inertia_apply is never echoed" and "derived, not invented" invariants
executable checks.

**Risk.** Low. The drop-list is the sensitive part and is version-sensitive;
the validation mode bounds the risk by failing loudly when the derivation no
longer reproduces the known pairs.

**Effort.** **medium**.

---

## f. Architectural debt found

The cycle surfaced several pieces of debt. Each is a small, independent cleanup;
none blocks the current results.

### f1. Three replicated quadrature defaults

**Problem.** The element-integral quadrature default is replicated in three
places: the element default, `msf_element_rule` (`calcul_force.cc`) and
`mc_element_rule` (`delete.cc`). They were switched to Gauss together in
`06e175d`, but a future change to the default must be kept in sync manually.

**Proposed change.** Resolve the rule from a single helper (e.g.
`element_quadrature_rule()`), with the two call sites reading it.

**Benefit.** One place to change; no silent divergence between the element
integral and the `msf`/`mc` paths.

**Risk.** Low (a mechanical refactor; A/B the corpora).

**Effort.** **cheap**.

### f2. Non-principal stress dofs (the staggered workaround)

**Problem.** `materi_stress` dofs are non-principal (`input.cc` sets no
`dof_principal`), so they are updated after the solve from the pre-solve
constitutive stress. The one-iteration lag was closed by a post-solve pass
(`dae06a3`), which is a workaround, not the structural fix. It also carries the
performance cost of §a.

**Proposed change.** Revisit the monolithic saddle-point mixed solve
(`mpc5_staggered_stress_lag.md` candidate (a)) **after** making the tie scheme
energy-consistent for shared nodes (the reason candidate (a) was refuted). If
that is not pursued, at minimum keep the post-solve pass as the documented
mechanism and apply §a.

**Benefit.** A structurally correct stress dof, no reporting lag, and the
removal of the double element pass.

**Risk.** High (corpus-wide refactor; the shared-node consistency must be solved
first).

**Effort.** **large**.

### f3. Reset scoping legacy branch

**Problem.** `control_reset_dof` keeps a legacy every-step behaviour for resets
whose index also carries a `control_timestep` (`control_reset_dof_scoping.md`).
The Professional rejects such inputs, so the branch is unreachable for
Professional-aligned models, but it is a GNU-only semantic that can surprise a
user.

**Proposed change.** Either reject reset+timestep at the same index (matching
the Professional) or document the branch prominently as a GNU extension and add
a warning.

**Benefit.** One clear reset semantic; no silent GNU/Pro divergence.

**Risk.** Low; the legacy branch is not exercised by any corpus test (A/B
byte-identical), but a legacy GNU input could rely on it.

**Effort.** **cheap**.

### f4. `ELEMENT_MATRIX_DELETE` dead record

**Problem.** `ELEMENT_MATRIX_DELETE` and its allocation in `top.cc`/
`database.cc` remain registered but are now unused (the Professional has no
equivalent; `delete_factor_element_matrix_scaling.md`). The deleted-element
path no longer stores a matrix snapshot.

**Proposed change.** Retire the record and its allocation, and the enum slot if
a renumbering of the database enum is acceptable (the headers
`tochnog.h`/`tochnog-mod.h` must stay in sync).

**Benefit.** Less dead code; no risk of a future contributor re-introducing the
diagonal-snapshot behaviour.

**Risk.** Low; the enum renumbering must be done carefully (the database enum is
shared between the two headers).

**Effort.** **cheap**.

### f5. `refine_l.cc` unreachable BAR3 branch

**Problem.** A second copy-paste defect in `refine_locally()`: the side-node
list inside the `for (iside=...)` loop has `else if (name==-BAR2)` where
`-BAR3` was intended (`refine_l.cc:249`), so the branch is unreachable and the
`-BAR3` local-refinement path would use a stale `ilist`. It is not exercised by
any corpus test (the only `control_mesh_refine_locally` bar model is `-bar2`).

**Proposed change.** Change the branch to `-BAR3` and add a focused `-bar3`
local-refinement test.

**Benefit.** The `-BAR3` local-refinement path becomes correct and tested.

**Risk.** Low; no corpus test currently exercises it.

**Effort.** **cheap** (fix) / **cheap–medium** (test).

### f6. Band width derived from `NODE_NODE`

**Problem.** The Bi-CG direct-LU retry band width is derived from the
`NODE_NODE` connectivity with an identity ordering
(`(1+diff)*nprinc`). It is correct but conservative; a tighter exact bound could
be derived from the element matrix unknown pairs (`solver_bicg_band_width.md`
§Pending refactorings).

**Proposed change.** Derive the band from the actual element unknown pairs (or
add a debug-build bounds assertion in the band fill to catch any future
out-of-band write).

**Benefit.** A tighter band (less memory/CPU) and a guard against the silent
corruption class that `6eb399d` fixed.

**Risk.** Low. A wrong bound is dangerous, but the debug assertion makes it
fail loudly.

**Effort.** **cheap** (assertion) / **medium** (exact bound).

### f7. `control_eigen` is GNU-only

**Problem.** `control_eigen` / `node_eigen` have no Professional equivalent
(neither the 25-10-2023 nor the 02-08-2026 binary contains the strings), so the
eigen path has no convergence reference; the `ddc9a5b` fix restored the original
GNU semantics, not a Professional alignment.

**Proposed change.** Keep it documented as GNU-only and add an explicit note in
the manual-user page (there is no Professional behaviour to converge to).

**Benefit.** Prevents a future "convergence" work unit from chasing a
non-existent reference.

**Risk.** None.

**Effort.** **cheap**.

---

## Summary

| id | improvement | benefit | risk | effort |
|---|---|---|---|---|
| a | hoist the post-solve stress pass — **DONE** | halve added assembly cost | low–moderate | cheap–medium |
| b | preconditioner; revisit the mixed solve | robustness; structural fix | medium–high | medium/large |
| c | clean the `build_safe.sh` full-mode file-checks | a clean CI gate | moderate | medium |
| d | reproducible corpora + stable runner + recalibration tooling | reproducibility, CI | low–moderate | medium |
| e | promote the `.dbs` reconstruction tooling | re-runnable, self-checking recovery | low | medium |
| f1 | single quadrature-rule helper | one place to change | low | cheap |
| f2 | principal stress dofs (mixed solve) | structural correctness | high | large |
| f3 | reset scoping | one clear semantic | low | cheap |
| f4 | retire `ELEMENT_MATRIX_DELETE` | less dead code | low | cheap |
| f5 | fix the `-BAR3` refine branch | correct, tested path | low | cheap |
| f6 | exact band width / bounds assertion | tighter band; loud failure | low | cheap/medium |
| f7 | document `control_eigen` as GNU-only | avoid a phantom work unit | none | cheap |

Highest value for the least risk: **c** (a clean build gate), **d** (a
reproducible corpus), and **a** (the performance cost). The principled
long-term fix is **b/f2**, but it depends on the shared-node tie consistency
and is a corpus-wide refactor.
