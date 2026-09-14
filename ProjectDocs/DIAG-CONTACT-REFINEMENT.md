# DIAG: contact/refinement corpus subgroup (examp6, examp17)

Date: 2026-09-14
Branch: `dev/contact-refinement` (from `6eb399d`)
Baseline commit: `6eb399d` (`dev/masin-kernel-surface`)
Corpus baseline: **187 PASS / 10 RUNFAIL / 0 PARSE** (197 tests)
Evidence: `/tmp/opencode/contact-20260914-184827/evidence/`
Scratch worktree: `/tmp/opencode/contact-20260914-184827/wt`

Binaries compared (all runs `LC_ALL=C`):

| id | binary | what it is |
|----|--------|------------|
| GNU-current | `.../contact-20260914-184827/wt/build/tochnog` | clean build of `6eb399d` |
| GNU-2014 | `/tmp/opencode/tn2014/build/tochnog` | F. Lorenzo, Sources.zip jan-2014 + our `hypo.c` |
| GNU-2014src | `.../contact-20260914-184827/orig/build/tochnog` | repo `73f0e70` (2026-08-02: 2014 sources + first modern-g++ build fixes); numerically identical to GNU-2014 |
| GNU-pre | `/tmp/opencode/oldrepo/build/tochnog` | 2026-09-02, before the `inertia_apply` alignment |
| Pro | `/tmp/opencode/tn_prof/tochnog_version_02-08-2026/src/tochnog` | Professional build 02-08-2026 |

## Result

| test | target | GNU-current | GNU-2014 | Professional | root cause | class |
|------|--------|-------------|----------|--------------|------------|-------|
| examp6 | node 259 x = 0.51 ±0.02 | 0.609375 (299 nodes) | **0.5078125 (274 nodes, rc=0)** | no own copy | regression from the `ELEMENT_DOF` stress/history slot fix `8f16dc0`, compounded by the `inertia_apply` default `0560dcf`; target tied to the residue-driven local-refinement node numbering | OUT-OF-SCOPE |
| examp17 | sigyy = -20506 ±100 | -21103.3 | -21330.5 | `validation_10`: -21031.0 (no p) / -20699.7 (with p); own target -20538.3 ±1e4 | stale target: no GNU or Professional binary reaches -20506 ±100 | OUT-OF-SCOPE |

Both are **OUT-OF-SCOPE** for a minimal source fix. The prior triage's
"pre-existing" label is correct for examp17 but **wrong for examp6**: examp6
PASSES on the 2014 baseline and on the repo's first modern build, so it is a
regression introduced by later project commits.

## 1. examp6 — regression from the `ELEMENT_DOF` slot fix + the inertia default

Input: unit square, one `-quad4`, von Mises plasticity with linear hardening,
`contact_geometry 0 -geometry_point 1`, two global `-h_refinement`, one
`control_mesh_split -yes`, `control_mesh_refine_locally 150 20.`, a restart and
three repeats. Target: the x coordinate of the remesh-created node 259.

### Evidence

- The sfnet corpus ships a **reference database** `examp6.dbs` (274 nodes, 508
  `-tria3`). Its node 259 is `(0.5078125, 0.9296875)`, inside the target band.
- GNU-2014 and GNU-2014src reproduce that database exactly: rc=0, 274 nodes,
  node 259 x = 0.5078125. So the target is valid for the project's source
  baseline.
- GNU-current: rc=1, 299 nodes, node 259 x = 0.609375.
- The target coordinate `(0.5078125, 0.9296875)` still exists in the
  GNU-current mesh, but at **node 274** (default) / **node 262**
  (`inertia_apply -yes`). The node is placed; only the numbering shifted.
- With `inertia_apply -yes` (the pre-`0560dcf` default): rc=1, 279 nodes,
  node 259 x = 0.53125 — misses the band by 0.00125.
- GNU-pre (2026-09-02) already gives 279 / 0.53125, i.e. it already carries the
  first regression.
- `git bisect` (clean build + run at each step): first bad commit **`8f16dc0`**
  ("feat(solver): estado sigma en equilibrio post-solve (B^T sigma = P)",
  2026-08-29); its parent `03cad89` passes (274 / 0.5078125).
- `8f16dc0` changes the `ELEMENT_DOF` slot layout for `nder>1`: examp6 uses
  `derivatives` (2D ⇒ `nder=4`), so the remembered integration-point
  stress/history ranges (`elem.cc`) and the constitutive write index
  (`materi.cc`) change. That changes the mapped material state across the
  local refinement/restart, hence the equation residue, hence which elements
  the residue-driven `control_mesh_refine_locally` selects.
- Reverting `8f16dc0`'s `materi.cc` + `elem.cc` hunks in the `6eb399d` tree does
  **not** restore the 2014 mesh (282 nodes / 0.546875): later commits change the
  same path again. A clean restore would need reverting the whole chain of
  deliberate solver/`ELEMENT_DOF` improvements.
- With `control_mesh_refine_locally` removed, GNU-2014src and GNU-pre produce the
  **same** mesh (41 nodes / 64 elements). The divergence is therefore exactly the
  residue-driven selection, not the global refinement/split.

### Why OUT-OF-SCOPE

`8f16dc0` is a correct fix: before it, the constitutive stress was written to the
unknown-number slot (`stres_indx/nder + j`) instead of the value slot
(`stres_indx + j*nder`), so for `nder>1` the element integration-point stress
stayed zero. Reverting it — or the later solver improvements — would regress the
`gforce7`/`gforce10`/`gforce13`/`qsri3d` family. The examp6 target is tied to the
pre-fix node numbering of a residue-driven adaptive remesh, which is inherently
solver-sensitive; there is no minimal source fix that satisfies it without those
regressions.

## 2. examp17 — stale target (the Professional's own copy uses another model and a 100x looser target)

Input: quarter of a circular sector (one `-quad4`, node 3 on
`geometry_circle 0`), elasticity, `contact_geometry 1 -geometry_line 2`, three
global `-h_refinement` plus one `-p_refinement`, `-updated_without_rotation`.
Target: `post_point_dof 0 -sigyy = -20506 ±100` at t=1.

### Evidence

- GNU-current: rc=1, -21103.3. Time series (dt=0.1, 10 steps): -3385.29,
  -6439.15, -9519.68, -12015.3, -14056.8, -16109.2, -17660.5, -18811.6,
  -19959.2, -21103.3.
- GNU-2014 and GNU-2014src (with the legacy `number_of_integration_points 9`
  record the old input needs): rc=1, -21330.5.
- GNU-pre (2026-09-02): rc=1, -21330.5.
- The **Professional cannot run the corpus file**: it does not know
  `-updated_without_rotation`, `bounda_unknown` or `geometry_set` (GNU dialect).
  The Professional's own equivalent is `test/validation/validation_10.dat`:
  - identical nodes/element/geometry/contact line, but `materi_displacement` +
    `group_materi_memory 0 -total_linear`, `bounda_dof`,
    `contact_target_geometry`, `-size_dev`, **no p-refinement**, and target
    `-20538.3 ±1e4` (100x looser);
  - Professional run: final `post_point_dof 0 = -21031.0`, rc=0;
  - Professional run with the corpus p-refinement (step 39) added: final
    -20699.7, rc=0.
- GNU on the same `-total_linear` + `materi_displacement` model: corpus + p =
  -20796.1; corpus without p = -21232.4. GNU-2014src + total_linear + p =
  -20946.7.
- The initial (`node_start_refined`) mesh is **byte-identical** between
  GNU-2014src and GNU-current (289 nodes, 64 elements): the divergence is the
  solution (formulation/solver), not the h/p refinement.
- Therefore the stored `-20506 ±100` is met by no available binary: the closest
  are the Professional's -20699.7 (194 outside the band) and the GNU
  `-total_linear` -20796.1 (290 outside). The target predates the project's 2014
  baseline (the 2014 GNU already gives -21330.5).

### Why OUT-OF-SCOPE

The target is stale relative to both the GNU 2014 baseline and the current
Professional; there is no source change that reaches -20506 ±100. Aligning the
input with the Professional's `validation_10` (the manual dialect) would change
the test and still land at ~-20700/-21030, not -20506. (The corpus file is
identical to the original `test/examp17.dat` shipped in the GNU tarballs.)

### Side finding (not fixed, out of scope)

The `contact_target_geometry` alias added in `contact.cc` only copies index 0
(`db_active_index(CONTACT_TARGET_GEOMETRY, 0, ...)`). The Professional's own
tests use index 1 (`contact_target_geometry 1 -geometry_line 2`), so the GNU
silently drops the contact and returns ~0 for `validation_10`. The corpus tests
all use the GNU name `contact_geometry`, so they are unaffected. Fixing it
properly requires moving the alias out of the parallel region
(`parallel_contact` runs through `parallel_sys_routine`; the `PUT` allocates —
see `manual-developer/contact_heat_generation.md`), so it is not a one-line
change.

## Validation

- Full 197 corpus with the `6eb399d` worktree binary (`final-corpus-run.sh`,
  `LC_ALL=C`): **187 PASS / 10 RUNFAIL / 0 PARSE** — identical to the stated
  baseline. No source change was made.
  - RUNFAIL: `examp6`, `examp14`, `examp17`, `examp24`, `ho_mech2`, `ho_mech3`,
    `ho_mech4`, `petsc1`, `pruebavel`, `taylor3`.
- Internal suite `scripts/build_safe.sh`: **16/16 runs OK** (`hypo1-4`,
  `gforce7`, `tsup_*`).
- No source files were modified.

## Conclusion

No minimal, source-grounded fix exists for this subgroup:

- **examp6** is a regression, but from a chain of deliberate, correct
  solver/`ELEMENT_DOF` improvements (first `8f16dc0`); the target is tied to the
  pre-fix residue-driven mesh numbering.
- **examp17**'s stored target is stale: it is not reproduced by the 2014
  baseline, the current GNU, or the Professional's own equivalent test.

Future work: examp6 needs either a target recalibration to the corrected remesh
numbering or a dedicated adaptive-remesh work unit; examp17 needs its
input/target aligned with the Professional `validation_10` semantics (a separate
convergence work unit), not a contact/refinement bug fix.
