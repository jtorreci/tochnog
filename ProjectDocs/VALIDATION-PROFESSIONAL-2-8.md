# Professional validation 2 and 8 — closure

Date: 2026-09-15
Branch: `dev/validation-2-8` (worktree `/tmp/opencode/val28-20260915-202314/wt`)
Baseline: `a41cda4` (`dev/masin-kernel-surface`), clean `build_safe.sh --clean`
(`TN_SKIP_SUITE=1`).
Evidence: `/tmp/opencode/val28-20260915-202314/`.

Binaries:

| binary | path | notes |
|---|---|---|
| GNU (this branch) | `wt/build/tochnog` | clean build from `a41cda4`, md5 `5d562a0b8fb8a9b0de9e79a608a36ece` |
| Tochnog Professional | `/tmp/opencode/tn_prof/tochnog_version_02-08-2026/src/tochnog` | static, build 02-08-2026 |
| GNU jan-2014 (F. Lorenzo) | `/tmp/opencode/tn2014/build/tochnog` | reference fork |
| GNU pre-merge (`oldrepo`) | `/tmp/opencode/oldrepo/build/tochnog` | earlier tree, sanity cross-check |

All runs with `LC_ALL=C`,
`LD_LIBRARY_PATH=<repo>/external-downloads/numlib-runtime`.

## Where the tests live

The `validation_*` tests are **not** in the sfnet corpus
(`external-downloads/sfnet/extracted/test`, 197 `.dat`, no `validation_*`)
and **not** in `validation-suite/test-2014/` (233 reconstructed tests, no
`validation_*`). They are the Professional's own validation suite, shipped
in `tochnog_version_02-08-2026/test/validation/validation_1..14.dat`
(extracted working copy:
`/tmp/masin-kernel-corpus-full.8pjO3w/input/validation/`). Both inputs and
the corpus are **untracked** (gitignored); this document records the exact
evidence.

## Result

| test | class | rc | value / check |
|---|---|---|---|
| `validation_2` | **FIXABLE (source, already fixed by `6eb399d`)** | 0 | `post_point_dof 1 -pres = -0.634519741323` = Pro `-0.6345197413352` (10 digits); target `-0.63 +/- 0.05` |
| `validation_8` | **FIXABLE (source, already fixed by `f244ac5`)** | 0 | 50/50 `mpc_node_number` + 8/8 `mpc_node_factor` records byte-identical to the Pro; no target record in the input |

Neither test needed a new source change: both were already closed by
commits that are ancestors of the `a41cda4` baseline. The prior notes
(`SEGUIMIENTO` lines 60/833) predate those commits and are **stale**.

## validation_2 — coupled 1D consolidation

Input: 1D bar, `materi_velocity`/`materi_displacement`/`materi_stress`/
`groundflow_pressure`, `inertia_apply -yes` +
`groundflow_consolidation_apply -yes`, 6 global h-refinements,
`control_timestep 20 delta_time 2.2`, target
`post_point_dof 1 -pres = -0.63 +/- 0.05`.

Root cause of the historical RUNFAIL: the Bi-CG direct-LU retry band
width in `so.cc solve()` was sized `diff*nprinc` instead of
`(1+diff)*nprinc`, so the pressure<->velocity consolidation coupling
(`nprinc=2`, dof distance 3 > band 2) wrote outside `mat`/`ldmat` and
corrupted the `dgbsv` retry. The default solver breaks the Bi-CG
(ill-conditioned consolidation system) and takes that retry, so the
result was wrong (`-0.7206` in `examp2`, the sfnet copy of the same
physics; `-0.688985` in the earlier `validation_2` note). With a forced
direct solver the assembled system was already correct
(`-0.634519741323`).

Fix: commit `6eb399d` (`fix(solver): correct Bi-CG direct-LU retry band
width`), ancestor of `a41cda4`. See
`manual-developer/solver_bicg_band_width.md` and
`CORPUS-PARSE-INIT-EXAMP2.md`.

Measured values (this work unit):

| binary | `post_point_dof 1 -pres` | rc |
|---|---:|---:|
| GNU `a41cda4` (default solver) | `-0.634519741323` | 0 |
| GNU `a41cda4` (`control_solver 0 -matrix_lapack`) | `-0.634519741323` | 0 |
| Professional 02-08-2026 | `-0.6345197413352` | 0 |
| GNU pre-merge (`oldrepo`) | `-0.634519842135` | 0 |
| GNU jan-2014 | n/a — rejects `inertia_apply` | 1 |

Conclusion: the stored target `-0.63` is **not** stale — the model gives
`-0.6345`, which is inside `+/- 0.05`. The honest reference value is the
Professional's `-0.6345197413352`. The `-0.688985` in the prior note was
the corrupt-retry artifact, refuted by both the current and the pre-merge
GNU.

## validation_8 — continuous metal forming (`-veln`)

Input: 2D eulerian (`mesh -fixed_in_space`, `convection_apply -yes`),
von Mises plasticity with hardening, boundary
`bounda_dof 2 -geometry_set 1 -veln` (zero normal velocity on the upper
walls). No `target_item`/`target_value`: the pass criterion is `rc=0`
(confirmed in the source: `miscel.cc` exits only when a target exists and
is missed).

Root cause of the historical PARSE stop: `-veln` was not implemented.
Manual Professional 6.22 (user manual, "bounda_dof") says the normal
direction comes from the geometrical entity (or `bounda_normal` for node
ranges) and that "internally in Tochnog a multi-point-constraint will be
generated". The manual does not specify which velocity axis becomes the
slave and which the masters — that is an implementation detail that was
deduced by measurement against the Professional and implemented in commit
`f244ac5` (`feat(bounda): -veln normal-velocity bounda_dof with generated
mpc ties`), ancestor of `a41cda4`:

- normal = `geometry()` first-match entity in record order (a
  `geometry_set` is scanned in order; corner nodes resolve to the first
  entity);
- one `mpc_node_number` per boundary node (ascending node order): slave =
  first velocity axis with a non-zero normal component, masters = the
  remaining axes, factor `-n_master/n_slave`; zero-component masters are
  omitted (an axis-aligned wall collapses to a masterless bound-to-zero
  record);
- the generated records carry `mpc_from_bounda -yes` and are regenerated
  only when the mesh fingerprint changes.

Fresh evidence (this work unit): the GNU `a41cda4` run reproduces the
Professional's records exactly — 50/50 `mpc_node_number` and 8/8
`mpc_node_factor` (values `-1.0e1`; the Pro prints scientific notation,
the GNU decimal — numerically identical). The `-veln` axis-choice
algorithm is therefore **implementable and implemented**; the prior
"not deducible without more data" note is stale.

Measured: GNU `rc=0`, wall `163.5 s`; Pro `rc=0`, wall `200 s`. Both
codes run the same 2000 steps over 192 elements, so the runtime is
inherent to the test, not a GNU inefficiency (a 20-step variant is
`1.6 s`, i.e. linear).

Field agreement (separate from the tie structure): the velocity field
matches the Professional within ~0.4% (max `|dvx,dvy| = 0.0042` over the
225 nodes), consistent with the `-veln` constraint being imposed
correctly. The post-solution stress/plastic-strain field diverges
(`post_point_dof_calcul` deviatoric stress `4.99` vs Pro `22.3`) — that is
the eulerian convection / plasticity family documented for `validation_1`
(the `convection_apply` term of the mixed velocity/stress assembly), a
separate work unit, not the `-veln` machinery. The test has no target, so
this does not affect its `rc=0` classification.

### Residual — the corpus timeout (not a source bug)

The test does not fit the 45 s corpus slot (163-200 s in both binaries),
so a corpus runner classifies it as `DEFERRED`/timeout (`rc=124`). This
is a **harness budget** matter: no source fix can make an inherently
2000-step, 192-element eulerian plastic run ~4.5x faster when the
Professional is equally slow. It is recorded here as the honest
integration note; the test passes standalone (`rc=0`).

## Verification

| check | result |
|---|---|
| `validation_2` standalone | rc=0, `-0.634519741323` (Pro `-0.6345197413352`) |
| `validation_8` standalone | rc=0, 50/50 MPC records identical to Pro |
| sfnet 197 A/B | **187 PASS / 9 RUNFAIL / 1 PARSE**, per-test diff vs the `a41cda4` predecessor run = **0** |
| validation-suite 233 A/B | **201 PASS / 32 RUNFAIL / 0 PARSE**, per-test diff = **0** |
| internal suite | **16/16** runs OK |

Note on the sfnet totals: this run classifies `pruebavel` as `PARSE`
(`rc=1`, unknown-record error), giving `187/9/1`. The task brief quotes
`187/10/0`; the only difference is that single classifier bucket. The
PASS count (187) and the per-test A/B (0 changes) are identical to the
`a41cda4` predecessor run (`/tmp/opencode/mpc5del-20260915-180036/
corpus-sfnet-final/final.tsv`). `pruebavel` is a known dead/malformed
test (`CORPUS-PARSE-INIT-EXAMP2.md`).

## Reproducibility

- Worktree: `git worktree add /tmp/opencode/val28-20260915-202314/wt a41cda4`;
  `external-downloads` symlinked, `validation-suite/test-2014/*.dat`
  copied (gitignored inputs).
- Clean build: `TN_SKIP_SUITE=1 ./scripts/build_safe.sh --clean`
  (`build-clean.log`).
- Per-test runs: `runs/gnu2-fresh/`, `runs/pro2/`, `runs/gnu8/`,
  `runs/gnu8-rc/`, `runs/pro8/`, `runs/gnu2014-2/`, `runs/gnuold-2/`.
- Corpus runner: `evidence/run-corpus.sh`; outputs
  `evidence/corpus-sfnet/final.tsv`, `evidence/corpus-vs/final.tsv`.
- Manual extraction: `evidence/validation.txt` (Professional
  `doc/validation.pdf`, "Validation 8: Continuous metal forming").
