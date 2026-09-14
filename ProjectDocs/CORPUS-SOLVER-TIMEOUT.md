# Solver / timeout corpus subgroup (sfnet 197)

Date: 2026-09-14
Branch: `dev/solver-timeout-family`
Baseline: `b83d907` (branch `dev/masin-kernel-surface`), clean build.
Corpus: `external-downloads/sfnet/extracted/test` (197 `.dat`, **untracked** —
the exact edit is recorded below since git does not track it).
Evidence: `/tmp/opencode/solvertimeout-evidence-20260914-152324/`
GNU baseline binary: `/tmp/opencode/solvertimeout-base-20260914-152324/build/tochnog`
(clean `build_safe.sh` build from `b83d907`).
GNU fixed binary: `/tmp/opencode/solvertimeout-20260914-152324/build/tochnog`
(clean `build_safe.sh` build from `b83d907` + the `refine_l.cc` fix).
Reference GNU 2014: `/tmp/opencode/tn2014/build/tochnog` (F. Lorenzo jan-2014
sources, built with the local lapack/blas runtime).
Professional: `/tmp/opencode/tn_prof/tochnog_version_02-08-2026/src/tochnog`
(build 2-8-2026) and its own test copies under `.../test/other/`.
All runs with `LC_ALL=C`, corpus timeout 45 s as in the corpus runner.

## Result

| corpus | PASS | RUNFAIL | PARSE |
|---|---|---|---|
| baseline `b83d907` | 184 | 13 | 0 |
| after this work unit | **186** | **11** | **0** |

Per-test diff: **2 RUNFAIL -> PASS** (`incnav3`, `refine5`), 0 regressions.
Suite `build_safe.sh`: **16/16**.

## Triage

| test | symptom | root cause | class | evidence |
|------|---------|------------|-------|----------|
| `refine5` | abort rc=134 | missing `else` in the `refine_l.cc` element-shape dispatch: a `-bar2` element falls through to `assert(name==-TET10)` | **FIXABLE (source)** | 2014 GNU also aborts rc=134 at the same assert; after the fix rc=0, `post_point_dof 0 -temp = 1` (target `1` +/- `5e-2`) |
| `incnav3` | `system (numerically) singular or not positive definite` | no velocity BC + material with only `elasti_compressibility` (no shear) + **no inertia** -> the momentum operator is a pure volumetric operator, numerically singular. The `0560dcf` Professional alignment flipped `options_inertia` default `-yes` -> `-no`, so this legacy model lost the mass term that regularised the operator | **FIXABLE (corpus input)** | 2014 GNU (inertia default `-yes`) passes with `node_dof 2 -tvely = -2.49461364553`; current GNU with `inertia_apply -yes` gives `-2.49461` (target `-2.5` +/- `1e-2`); the Professional's own `test/other/incnav3.dat` explicitly contains `inertia_apply -yes`; the Professional's PARDISO on the no-inertia GNU input fails with `Zero pivot` |
| `examp14` | timeout at 45 s | 4 global h-refinements + 250 steps; honest Bi-CG does **80** iterations (error `9.9e-17`) where the 2014 GNU stopped at **8** (stagnation exit, error `9.4e-13`). Runtime 49 s (GNU) vs 19 s (2014). The stored target `15.131` is **stale**: no GNU build reaches it (current `9.71888`, 2014 `9.61289`) | **OUT-OF-SCOPE** | see "Timeouts" below; no Professional copy of `examp14` exists, and the GNU `options_mesh`/`bounda_unknown` dialect is rejected by the Professional binary |
| `examp24` | timeout at 45 s | 150 `control_repeat` cycles of `control_mesh_new_mesh` (0.25) + one step. The cost is **mesh generation**, not the solver: the remesh alone (4860 tet4) takes **12.8 s** (GNU) / **7.1 s** (2014); the step is 37 ms. Full run: >300 s (GNU, killed) / >17 min extrapolated (2014). No target record in the input | **OUT-OF-SCOPE** | 150 x ~13 s ~ 32 min GNU; the 45 s corpus budget is unrealistic for this test; 2014 is ~1.8x faster but still ~18 min |
| `ho_mech3` | timeout at 45 s | 40 `control_repeat` cycles of a macro mesh + 100 steps. Both binaries **diverge** (post_node_result ~1e15 vs target 22000) | **OUT-OF-SCOPE** | GNU 302 s rc=1 (`post_node_result 1 = -7.19e14`), 2014 164 s rc=1 (`-1.48e15`); the divergence is present in the reference GNU, so it is a model/numerics limit, not a regression |
| `petsc1` | `Error: set PETSC_USE to 1 in tnpetsc.h` | PETSc is not a build option in this environment | **OUT-OF-SCOPE (infrastructure)** | `tnpetsc.h` has `PETSC_USE 0`; no `libpetsc`/`petsc.h` in the environment; `so_petsc.c` is entirely inside `#if PETSC_USE`; the Professional 02-08-2026 source tree ships no `tnpetsc.h` either |

## 1. `incnav3` — ninth victim of the `inertia_apply` default flip

`incnav3` is the nearly-incompressible Navier-Stokes example: one `quad4`,
all four nodes with prescribed temperature 0.5, `force_gravity 0. 1.`, and a
material with only `elasti_compressibility 1.0` (no `elasti_young`, no
Poisson -> no shear stiffness). There is **no velocity boundary condition**.

The momentum operator is therefore the volumetric operator
`K_vol = integral (div N_i)(div N_j)` alone, which is numerically singular
(its null space contains every divergence-free velocity field plus the rigid
body modes). The old GNU ran it because the historic `options_inertia`
default was `-yes`: the mass term `M/dt^2` adds to the diagonal and makes the
operator positive definite. Commit `0560dcf` ("align inertia_apply and
consolidation defaults with the Professional", an ancestor of `b83d907`)
flipped the default to `-no` (Professional manual 6.785), so this legacy
model lost the regularising mass term.

Measured `node_dof 2 -tvely` (target `-2.5` +/- `1e-2`):

| run | value | rc |
|-----|-------|----|
| 2014 GNU (inertia default `-yes`) | -2.49461364553 | 0 |
| current GNU, no inertia (corpus file) | Bi-CG breakdown, then LU retry fails | 1 |
| current GNU + `inertia_apply -yes` | **-2.49461** | 0 |
| Professional 02-08-2026 on its own `test/other/incnav3.dat` (has `inertia_apply -yes`) | -4.990415335463 (its own target `-5.0`) | 0 |
| Professional 02-08-2026 on the no-inertia GNU input (dialect translated) | `Pardiso internal error number -4 / Zero pivot` | 1 |

The Professional's own copy adds `inertia_apply -yes` and also changes the
target to `-5.0` and the memory to `-updated_linear` (a different physics
variant: it lets `rho = rho_0 (1 - expansion T)` change the gravity term).
The GNU variant keeps its own target `-2.5`; with `inertia_apply -yes` the GNU
reproduces its own 2014 value `-2.49461`, which is inside the stored
tolerance. This is the same class as the eight stale inputs fixed in
`CORPUS-STALE-ALIGNMENT.md` section 1; `incnav3` was missed there because its
symptom is a singular-solver abort rather than a wrong steady value.

**Fix applied to the corpus input** (the file is gitignored): insert

```
inertia_apply -yes
```

immediately after the `end_initia` line of `incnav3.dat`. No source change.

## 2. Timeouts — where the time goes

### `examp14`

Input: `control_mesh_refine_globally` at indices 2, 3, 4, 5 (four h-refinements)
followed by `control_timestep 20 2.e-2 5.` (250 steps).

| binary | wall | `time_calculation` | `global_solver_iterations` | `global_solver_error` | target value |
|--------|------|--------------------|----------------------------|------------------------|--------------|
| current GNU | 48.95 s | 47 s | 80 | 9.90e-17 | 9.71887788387 |
| 2014 GNU | 19.14 s | 19 s | 8 | 9.39e-13 | 9.61289226266 |

The target stored in the corpus is `15.131`; **neither** GNU build reaches it
(difference > 55%), so even a sub-45 s runtime would still be a RUNFAIL. The
2.5x runtime difference is the honest Bi-CG converging to `1e-14` relative
(the 2014 GNU stopped at 8 iterations through the removed stagnation exit).
The Professional has no copy of `examp14` and rejects the GNU `options_mesh`
dialect, so its runtime is not measurable. Classified OUT-OF-SCOPE: the
timeout masks a stale, unreachable target, and the extra iterations are the
deliberate honest-solver behaviour, not a defect.

### `examp24`

Input: `control_mesh_new_mesh 10 0.25`, `control_timestep 100 1.e-6 1.e-6`,
`control_repeat 200 150 10` -> **150** remesh + step cycles. No target record
(a `rc=0` finish is the pass condition).

| run | wall | elements |
|-----|------|----------|
| mesh generation only (0.25) | 12.8 s (GNU) / 7.1 s (2014) | 4860 |
| one full cycle (`control_repeat ... 1 10`) | 18.9 s (GNU) | 4860 |
| mesh size 0.5 | 1.0 s | 900 |
| mesh size 1.0 | 0.25 s | 324 |
| full 150 cycles | >300 s (GNU, killed at 300 s) / ~18 min (2014 extrapolated) | 4860 |

The cost is dominated by `control_mesh_new_mesh` (mesh generation): the remesh
alone is 12.8 s while a full remesh + step cycle is ~14 s (marginal over the
first cycle), so the solve is ~1.5 s (`global_solver_iterations 10`).
The current GNU is ~1.8x slower at remeshing than the 2014 reference, but
neither is anywhere near the 45 s corpus budget: the test is intrinsically
long. Classified OUT-OF-SCOPE.

### `ho_mech3`

Input: `control_mesh_macro`, `control_mesh_new_mesh 20 8.`,
`control_timestep 40 1.e-2 1.` (100 steps), `control_repeat 70 40 20`
(40 cycles).

| binary | wall | target value (`post_node_result 1 -vely`, target 22000) |
|--------|------|--------------------------------------------------------|
| current GNU | 302 s | -7.19e14 |
| 2014 GNU | 164 s | -1.48e15 |

Both builds diverge (values ~1e14-1e15), so the stored target is not
reachable. Classified OUT-OF-SCOPE (model/numerics, present in the reference
GNU; the current GNU is ~1.8x slower).

## 3. `petsc1` — infrastructure

`tnpetsc.h` defines `PETSC_USE 0`; `so_petsc.c` is entirely inside
`#if PETSC_USE`; the error `Error: set PETSC_USE to 1 in tnpetsc.h` is the
deliberate stub. No PETSc headers or libraries are installed in this
environment (`ldconfig` finds no `libpetsc`, no `/usr/include/petsc*`), and
the Professional 02-08-2026 source tree ships no `tnpetsc.h` either.
Classified OUT-OF-SCOPE (build infrastructure), not a source defect.

## 4. Exact corpus edit (untracked file)

`incnav3.dat` — one line inserted immediately after the `end_initia` line:

```
+ inertia_apply -yes
```

No other corpus file is touched. The original is preserved in the evidence
dir under `originals/`.

## 5. Reproducibility

- Corpus runner: `/tmp/opencode/corpus-run.sh <binary> <corpusdir> <out> <work> 45`.
- Baseline: `b83d907` binary + original corpus.
- Fixed: `b83d907` + `refine_l.cc` fix + corpus copy with the `incnav3` edit.
- Per-test classification files: `corpus-baseline2.txt`, `corpus-fixed.txt`.
- Suite: `scripts/build_safe.sh` clean build + internal suite, 16/16.
