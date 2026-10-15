# solve_iterative_bicg — honest stopping criteria + CG for the SPD system

**File**: `so_bicg.cc` — the default iterative linear solver of the GNU
(`control_options_solver -matrix_iterative_bicg`, the default per timestep
in `top.cc`).

**Lot**: A+B of `DIAG-SOLVE-MIXTO.md` (2026-08-28). Fixes the dishonest
stopping criteria (A) and adds plain CG for the symmetric system (B).

## What changed

### 1. Honest stopping criteria (A)

The old loop had four exits, three of which reported *success* without
convergence:

| Old exit | Old behavior | New behavior |
|---|---|---|
| `|dAd| < 1e-16` (breakdown) | `ready = 1` → success, possibly with `x = 0` | honest failure (singular / not positive definite), unless the residual is already within tolerance (lucky breakdown) |
| `|r1·r2| < 1e-16` (bi-Lanczos collapse) | `ready = 1` → success with a partial solution | honest failure (non-converged residual) or convergence (residual within tolerance) |
| `|error − last_error| < 0.1·check_error` (stagnation) | `ready = 1` → success with `error ≫ check_error` (measured: 4.7e-7 vs 1e-12) | **removed** — CG/Bi-CG residuals are non-monotone; a flat step is normal and the iteration continues |
| `iter == max_iter` | failure (with `control_solver_bicg_stop` semantics) | unchanged, but now reports the honest relative residual |

**Convergence test** (the only success exit): the residual of the current
iterate, recomputed from the iterate on every pass,

```
error = |b_hat - A_hat x_hat|^2        (preconditioned system)
converged  <=>  error < check_error
check_error = max( bicg_error * |r0|^2, bicg_error_minimum )
```

with `bicg_error` default `1e-14` (`control_options_solver_bicg_error`,
overwritten by the global `solver_bicg_error`) and `bicg_error_minimum`
default `1e-16` (`control_options_solver_bicg_error_minimum`); the loop also
stops when `error < eps_stop = EPS_TMP = 1e-16` (code defaults in
`solve_iterative_bicg`, `so_bicg.cc`; measured `check_error 4.42565e-15` for
`|r0|^2 = 0.442565` and `1e-16` for smaller `|r0|^2`; this page said
`1e-10` / `1e-12` before F4 T6f C74). This is a **relative** residual test
(`|r|/|r0| < sqrt(1e-14) = 1e-7` when the floor does not bind, with an
**absolute** floor `|r|^2 >= 1e-16` for tiny loads; the §21 eliminated tie
systems tighten both to `1e-24` / `1e-26`) — it is not the "error decrease"
heuristic, and it is evaluated on the residual of the solution that is
actually returned (the iterate is not updated after the check).

The failure messages print the initial error, the final error, the
**unpreconditioned relative residual** `|b − A x|/|b|` (computed by
`solve_iterative_bicg_real_residual`) and the iteration count.

`control_solver_bicg_stop -no` (manual Professional 6.376) now also covers
the breakdown exits (continue with a warning + the current solution),
consistent with its max_iter semantics.

Record lookup (so_bicg.cc, 2026-09-29): combined precedence, last read
wins:

1. `solver_bicg_stop` (global);
2. `control_solver_bicg_stop` with index **0**;
3. `control_solver_bicg_stop` with the control index of the running
   timestep (`ICONTROL`, skipped when it is 0).

Measured on the Professional 02-08-2026 with a singular bicg system
(elasti7 mesh without supports, `solver -matrix_iterative_bicg`,
timestep index 20 unless stated): the control record overrules the global
one (manual 6.376), but ONLY the index-0 control record is read, for
every timestep; a record with the timestep's own non-zero index is
ignored. The manual documents a per-timestep index, so the project owner
judged this a Professional bug (same class as the index-0
`control_materi_elasti_k0` guard fixed in the GNU, R8) and chose the
combined precedence (decision 2026-09-29): the index-0 record keeps its
Professional meaning (inputs that rely on `control_solver_bicg_stop 0
...` behave the same), and a record at the timestep's own index is
honoured on top of it. This is a deliberate divergence from the
Professional only when such a non-zero index record exists (BACKLOG
C13). "stops" = Professional "Error: the BICG solver did not converge" /
GNU "Error: the iterative solver ..."; "continues" = no error (Pro) /
"Warning: ... continuing" (GNU):

| control_solver_bicg_stop | solver_bicg_stop | Professional | GNU `47b6fda` (index 0 only) | GNU combined |
|---|---|---|---|---|
| 20 -no | — | stops | stops | **continues** (divergence) |
| 0 -no | — | continues | continues | continues |
| 10 -no | — | stops | stops | stops |
| 0 -yes, 20 -no | — | stops | stops | **continues** (divergence) |
| 0 -no, 20 -yes | — | continues | continues | **stops** (divergence) |
| 0 -no (timestep 0) | -yes | continues | continues | continues |
| 0 -yes (timestep 0) | -no | stops | stops | stops |
| 0 -no (timestep 0) | — | continues | continues | continues |
| — | — | stops | stops | stops |

Extra probes, identical on the three binaries: `0 -no, 20 -no` →
continues; global `-no` only (timestep 20) → continues; nothing at
timestep 0 → stops. The GNU before `47b6fda` (active-index check at 0,
then a read at the current control index, then the global record
overwriting it) matched none of the mixed rows. Scratch inputs only (the
GNU "stop" in this singular case does not abort the run — it prints the
error and continues with rc 0 — so no self-checking target distinguishes
both paths; see BACKLOG C12).

Regression of the combined precedence (clean build, vs the `47b6fda`
binary): build suite 16/16 runs OK with the same 70 file-check failures
by name; validation-suite 233 PASS / 1 RUNFAIL, sfnet 187/10, test-2026
23/23; every `.dbs` identical after masking timestamps and denormals
(`examp14` hit the 45 s timeout under load — 46 s when run alone, rc 1 as
before, `.dbs` identical). No corpus input has a non-zero-index
`control_solver_bicg_stop` record.

### 2. CG for the symmetric system (B)

`solve_iterative_bicg_symmetric()` measures, per solve, whether the
assembled system is symmetric: the global matrix is symmetric iff every
element matrix is symmetric (the `NODE_LHSIDE` diagonal block is diagonal;
the element matrices are stored dense, both triangles, in `elem.cc`). For
every element and every pair of free dofs `(i,j)` it computes
`|Aij − Aji|/max(|Aij|,|Aji|)` and classifies the system as symmetric when
the maximum over all pairs is below `EPS_SYMMETRY = 1e-8`.

- **symmetric** → plain **CG** (one direction vector; `Ad2` accumulation
  skipped — the matvec cost is halved). This is the natural method for the
  velocity-only SPD system (`materi_velocity`, `temp`, `pres`); on a
  symmetric matrix the old Bi-CG reduced to CG anyway.
- **non-symmetric** → **Bi-CG** with the honest criteria. Measured
  non-symmetric systems: the beam formulation (the translation↔rotation
  blocks carry different `dtime` factors — genuine non-symmetry of the
  assembled matrix), and the plastic-slip 3D interface system (measured
  `max_rel_asym` up to 1.0; `dᵀAd < 0` during the iteration → indefinite).

CG's breakdown exit (`|dAd| < 1e-16` with a non-converged residual) reports
a genuinely singular or non-positive-definite system — the old code
reported `x = 0` as success (measured on the CCW-mesh case:
`dAd = -8.7e-23` with the load in the null space).

### 3. Primal residual monitor (correctness fix inside (A))

The residual accumulated by `solve_iterative_bicg_sys`/`_element` is now
the **primal** residual `b_hat − A_hat·x_hat` (accumulated at the row
index with the column solution). The old accumulation measured the
**transpose** residual `b_hat − A_hatᵀ·x_hat`, which equals the primal one
for symmetric matrices but **never vanishes at the solution of a
non-symmetric system** — the beam solve could never converge by the error
test and only "finished" through the false breakdown exit. For symmetric
systems the two accumulations are identical, so the behavior of every
symmetric solve is unchanged (verified: suite 199/199, converged-path
values byte-identical).

### 4. Assembled CSR operator (F4 T9, 2026-10-10)

The products of `solve_iterative_bicg_sys()` run on a CSR matrix
assembled once per solve (element-by-element fallback with mpc ties):
see [solver_bicg_csr_operator](solver_bicg_csr_operator.md).

## Measured behavior (evidence)

- **3D cantilever hex8 ×1, tip load** (the `A·b ≈ 0` family): the old code
  stopped at iteration 1 (`error 4.72727e-07`, flat-residual identity,
  RC=0 without solving). The new CG iterates through the flat step and
  converges: `iterations 0..3`, `final error 4.2e-37`,
  `relative residual 9.4e-16` — the DIAG's predicted "keep iterating"
  behavior.
- **3D cantilever hex8/hex27 with real load + MSF force output**:
  both now solve (RC=0, CG path). hex8 ×8: clamp moment `0.26·P·L`
  (1-in-thickness shear lock, same family as the 2D `0.2315×`); hex27 ×8:
  `mom ≈ 1.1·P·(L−x)`, `she ≈ P` — load-based 3D validation is unblocked.
- **2D quad4 fixed point (the KEY question)**: sweeping the outer
  staggered iterations 1..32 with the honest solver gives the SAME fixed
  point as the DIAG: plain `0.0185184` (0.23148×), SRI `0.0185298`
  (0.23162×) — **the fixed point is a property of the staggered scheme,
  not an artifact of the dishonest solver** (the 2D inner solve was already
  converging honestly). Fix C/D is still required for the section-moment
  deficit.
- **Road model (concrete 1-in-thickness on soil, plane strain,
  fixed-fixed, distributed load)**: now runs RC=0 with converged solves,
  but the statics check still fails: `|M_end| + M_center = 0.001435` vs
  `pL²/8 = 0.01` → `0.143×` (scheme fixed point, locked).
- **Beam tests (bmom family)**: the beam system is non-symmetric; with the
  primal monitor the Bi-CG converges to the exact solution through the
  normal path (the old code only "passed" via the false breakdown exit).
- **iface_3d_slip**: the plastic-slip 3D interface system is non-symmetric
  and indefinite; the honest Bi-CG breaks down. The test opts into the
  documented `control_solver_bicg_stop -no` behavior (continue with a
  warning) — the test verifies interface physics, not solver convergence.
- **Timing**: no test approaches the 120 s `build_safe.sh` timeout (the
  slowest is `msf_sheet3d` at ~1.4 s).

### Bi-CG breakdown and the direct band LU retry (F4 T6f C74, 2026-10-13)

**Breakdown test** (`so_bicg.cc`, Bi-CG branch of `solve_iterative_bicg`):
`|r1·r2| < EPS_TMP` or `|d2ᵀ A d1| < EPS_dAd`, both `1e-16`, absolute.
With `solver_bicg_stop -no` the iteration stops with the current iterate
and the warning "Warning: the Bi-CG iteration broke down (r1.r2 or d^T A d
~ 0, continuing)" (unchanged). Otherwise `solve_iterative_bicg()` returns
0 and sets the global flag `solve_iterative_bicg_breakdown`, keeping the
initial error and the relative residual for the diagnostic.

**LU retry** (`so.cc`, `solve()`, after the `bicg_solver` call): a failed
Bi-CG solve is re-solved with the direct band LU (`solve_band_retry`,
LAPACK dgbsv in a Cuthill-McKee ordering, F4 T7).

- breakdown + LU success: ONE line per solve, `Warning: Bi-CG breakdown,
  system solved with the direct band LU`;
- breakdown + LU failure: the full diagnostic
  (`solve_iterative_bicg_breakdown_report()`: "Error: the Bi-CG iteration
  broke down. ... Probably you have a wrong input file ...") and then the
  LAPACK error, exit;
- the displacement-primary second solve (no LU retry): the full diagnostic
  before its "did not converge" warning;
- a max_iter failure keeps its own error block and the one-time "Warning:
  Bi-CG did not converge - retrying with direct LU".

Before C74 the full "Error ... wrong input file" block was printed inside
`solve_iterative_bicg()` for EVERY breakdown, also when the LU retry then
solved the system exactly: 599 times on `validation1_steady1`, 38 of 956
corpus runs (Professional `test/other`, sfnet, test-2026, test-2014).
Numerical results are unchanged (corpus A/B `.dbs` identical; test
`bicg_breakdown_lu_retry1`, target = the SuperLU solution).

**Why it breaks down** (validation_1, traced): near the steady state the
right-hand sides are tiny (`|r0|^2` 1e-12 .. 1e-8 for most solves), so
the absolute 1e-16 thresholds fire while the scale-relative quantities are
healthy: `|r1·r2|/(|r1||r2|)` ~ 1e-3 .. 1e-4 (min 1e-8),
`|d2ᵀAd1|/(|d2||Ad1|)` ~ 1e-2 .. 1e-4. These are false breakdowns; the
exact LU retry happens to give a MORE accurate solution than a converged
Bi-CG, whose absolute 1e-16 floor binds on tiny systems. A scale-relative
breakdown test alone (measured: no breakdown on validation_1, -10 % wall)
therefore loses accuracy (examp16 `node_dof` vs SuperLU 0.43 % -> 2.45 %;
validation_1 with the inflow scaled 1e-5: 8 %); it must come with a
relative convergence floor. Resolved by C77 (next section): the
breakdowns are gone, with the relative stop.

### Relative stop and guarded direct retry (F4 T6f C77, 2026-10-15)

Owner decision 2026-10-15. Implemented in `so_bicg.cc`
`solve_iterative_bicg()` and `so.cc` (`solve()`, `solve_superlu_retry()`,
`solve_band_retry()`, `solve_retry_visit()`, `solve_memory_room()`).

**Stop tests** (each loop pass, on the recomputed residual
`error = |P(b - A x)|^2`, x0 = 0 so `|r0|^2 = |Pb|^2`):

1. convergence: `error < bicg_error*|r0|^2`, or `error < bicg_error_minimum`
   ONLY when the user gives `options_solver_bicg_error_minimum` (global,
   read first) or `control_options_solver_bicg_error_minimum` (control
   index, wins); `eps_stop` and the floor 1e-26 stay for the §21 tightened
   tie systems (`mpc_elim_tighten`); `error == 0` (b = 0) converges;
2. round-off floor: `error < (BICG_FLOOR_FACTOR*DBL_EPSILON*|A|)^2 |x|^2`
   (iter > 0), `|A|` = the largest `|A d|/|d|` seen (`d1`, `Ad1`);
3. breakdown: Bi-CG `|r1·r2| <= BICG_BREAKDOWN_REL |r1||r2|` or
   `|d2·Ad1| <= BICG_BREAKDOWN_REL |d2||Ad1|`; CG `|dAd| < EPS_dAd_REL*error`
   (unchanged);
4. `iter == max_iter` (unchanged, own error block).

A floor or breakdown stop (`stop_kind` 1 / 2) is accepted when the true
unpreconditioned `|b-Ax|/|b| <= BICG_FLOOR_ACCEPT`; otherwise the routine
returns 0 with `solve_iterative_bicg_breakdown = 1` (C74 semantics:
`solver_bicg_stop -no` keeps the iterate with a warning). There is no
stagnation window and no best-iterate return: the prototype measured that
they cost accuracy (fallbacks 26 -> 10958 with a 10-iteration window).
The CG takes the same relative stop (measured safe, see below).

**Direct retry** (`so.cc`): SuperLU first (`solve_superlu_retry`): the
compressed-column matrix is built from the same element entries as the
iterative products (`solve_retry_visit`, mpc tie expansion included) plus
`NODE_LHSIDE`, duplicates merged per column, `|diag| < EPS_MAT -> EPS_MAT`
(as the band LU and the `-matrix_superlu` assembly); dgssv runs under the
`-matrix_pardiso` memory guard (`solve_pardiso_memory_guard`), and after a
failed factorization (it leaks its workspace) the run skips this retry
(`solve_superlu_retry_failed`). Then the band LU, refused before the
allocation when `8 ldmat n > BAND_RETRY_MEMORY_FRACTION * room`
(`solve_memory_room()` = min(MemAvailable, soft RLIMIT_AS - VmSize)):
"Error: the direct band LU retry does not fit in memory" (equations,
band, MB needed / available), return `SOLVE_BAND_RETRY_NO_MEMORY`, exit.
Messages: "Warning: iterative solver breakdown, system solved with the
direct SuperLU" / "... direct band LU" (one per solve), "Warning: the
iterative solver did not converge - retrying with a direct solver" (once);
`solve_iterative_bicg_breakdown_report()` says "Error: the iterative solver
broke down" (CG and Bi-CG).

**Measurements** (scratch instrumented builds, 86 corpus inputs with a
Bi-CG solve, reference = the same input with `solver -matrix_superlu`;
local, 2 jobs):

| variant | direct retries | iterations (Bi-CG / CG) | solver s | wall s | closer to SuperLU (10x) | farther |
|---|---|---|---|---|---|---|
| HEAD `6d3e19e` (absolute floors) | 7378 band LU | 1.34 M | 57.6 | 275 | - | - |
| relative stop (this fix) | 13 SuperLU (0 band) | 2.22 M / 15.7 k | 69.8 | 285 | 19 | 0 |

Large meshes (`large/gen.py`, 1 job, median wall): C23 (55296 dof, Bi-CG)
HEAD 91.9 s / 2.66 GB with 2 band-LU retries (73.7 s), 6.8e-7 off
SuperLU -> 18.6 s / 0.45 GB, 8.7e-12; C36 (202612 dof) HEAD bad_alloc ->
83.2 s / 1.8 GB; C13 4.70 -> 2.90 s. CG: E15 (12288 dof) 127 -> 295
iterations, 1.85 -> 1.99 s, 1.4e-6 -> 4.0e-11 off SuperLU; E25 (52728
dof) 9.84 -> 10.49 s, 9.3e-7 -> 1.6e-11. Corpus A/B on omen (965 runs vs
`c79_runs`): no classification change except `ground11_nonsaturated`
(timeout -> rc 1 = its HEAD rc when not timed out) and `truss13` (rc 1 in
both; HEAD singular system accepted on the absolute floor and printed no
`element_truss_force`, now the breakdown is reported and SuperLU / band
LU refuse the singular matrix); of the 84 runs whose `node_dof` moved
more than 1e-5, 63 moved toward SuperLU, 10 stayed at the same distance
(path-dependent plasticity, `tsup_*`, `examp16`), 1 moved away
(`dynamic6`, diverging in both, rc 1), 10 have a different refined mesh
(`refine*`, `delete4`, `examp7`, `examp15`, `ho_othr2`, `tslv_ver`; rc
unchanged).

**Non-symmetric algorithm** (owner item 2, measured; the decision is
pending, see BACKLOG C77): BiCGStab (two products, no `A^T`) and
GMRES(30) were measured against the relative-stop Bi-CG on the same 86
inputs: BiCGStab 71 SuperLU retries (29 solves at max_iter: `distri3`,
`beam2d_2`, `bimet1`), 2 inputs farther from SuperLU (`tquad6_iface`
5e-8 -> 7.5e-3, `ground19_water_under_dam` 6.5e-7 -> 7.6e-6), solver
68.3 s; large meshes 9-13 % less solver time (C36 29.0 -> 25.9 s, wall
83.2 -> 79.4 s, same solution to 8e-12). GMRES(30): 379 max_iter solves,
4.8x the products, `ground11_nonsaturated` timeout. The classical Bi-CG
stays until the owner decides.

## External features / libraries used

None new — plain C++ and the existing `parallel_sys_routine` machinery.
The symmetry check is one serial pass over the element matrices
(≈ one matvec), run once per solve.

## Hardcoded parameters

- `EPS_dAd_REL = 1e-14` (CG breakdown, relative to `|r|^2`); `EPS_TMP =
  1e-16` and `EPS_dAd` are no longer stop thresholds (C77).
- `EPS_P = 1e-10` (diagonal guard, unchanged).
- `EPS_SYMMETRY = 1e-8` (symmetry classification threshold; exact symmetric
  physics measures ~1e-16, genuine asymmetries measure 1e-3..1).
- `MAX_PUKNWN_SYMMETRY = 32` (skip the check — fall back to Bi-CG — for
  exotic dof counts; the dense per-element buffer is `(nnol·npuknwn)²`).
- `max_iter = 10·solve_nlocal` (unchanged).
- `check_error = bicg_error·|r0|²` (`bicg_error` default 1e-14), raised to
  `bicg_error_minimum` only when the user record exists (C77).
- `BICG_FLOOR_FACTOR = 100` (round-off floor `100 eps |A| |x|`),
  `BICG_BREAKDOWN_REL = 1e-14` (Bi-CG breakdown), `BICG_FLOOR_ACCEPT =
  1e-6` (true relative residual accepted at a floor / breakdown stop);
  prototype sweeps: breakdown 1e-12 and floor factor 10 give the same
  corpus result (24-26 retries).
- `so.cc`: `BAND_RETRY_MEMORY_FRACTION = 0.7` (band LU budget),
  `PARDISO_SUPERLU_MEMORY_FRACTION = 0.7` (SuperLU retry guard),
  `EPS_MAT = 1e-8` (direct diagonal regularization).

## Pending refactorings / known limitations

- The consolidation coupling block (materi_stress + groundflow_pressure,
  `materi.cc:619` vs `groundfl.cc:116`) has a sign mismatch between the
  off-diagonal blocks (measured `max_rel_asym` up to 1.0) — the coupled
  system is non-symmetric by construction and takes the Bi-CG path.
- `control_solver_bicg_stop`: combined precedence (global, index 0,
  timestep index); diverges from the Professional's index-0-only reading
  when a record with the timestep's non-zero index exists (see above,
  BACKLOG C13).
- The staggered fixed point (`0.2315×` 2D, `0.143×` road) is the scheme's;
  the monolithic mixed solve (fix C) or a regularized staggered update
  (fix D) is the definitive repair for the user case.

## Differences with the Professional version

- **Iterative solver stop** (F4 T6f C74, C77; DIVERGENCES C77, class A,
  owner decision 2026-10-15). GNU: relative stop, scale-relative
  breakdown test, no default absolute floor; no false breakdown on
  `validation_1`. Professional 02-08-2026: default PARDISO (direct); its
  `solver -matrix_iterative_bicg` stops on an absolute floor (1e-5 scaled
  `validation_1`: Pro Bi-CG velx 8.9888e-7, PARDISO 9.77518e-7, GNU
  9.775171e-7). MKL PARDISO is not integrated (SEGUIMIENTO Descarte).

- **`control_solver_bicg_stop` at the timestep's own index** (deliberate
  divergence, project decision 2026-09-29; BACKLOG C13; SEGUIMIENTO
  Descarte). GNU: Reads, last one wins: `solver_bicg_stop`, index-0 control
  record, control record at the running timestep index (manual 6.376).
  Professional 02-08-2026: Reads only the index-0 control record for every
  timestep; a record at the timestep's non-zero index is ignored. Status:
  Deliberate divergence, project owner; GNU follows the manual. Also
  Pro-suspect. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **Bi-CG breakdown on a singular system** (known difference, not decided;
  BACKLOG C12). GNU: Prints the breakdown error each iteration, finishes
  with rc 0 and a `.dbs`. Professional 02-08-2026: Prints "the BICG solver
  did not converge" and stops without a `.dbs`. Status: Open, LOW. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
