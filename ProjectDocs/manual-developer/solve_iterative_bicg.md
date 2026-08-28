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

with `bicg_error` default `1e-10` (`control_options_solver_bicg_error`) and
`bicg_error_minimum` default `1e-12`
(`control_options_solver_bicg_error_minimum`). This is a **relative**
residual test (`|r|/|r0| < sqrt(1e-10)` when the floor does not bind, with
an absolute floor for tiny loads) — it is not the "error decrease"
heuristic, and it is evaluated on the residual of the solution that is
actually returned (the iterate is not updated after the check).

The failure messages print the initial error, the final error, the
**unpreconditioned relative residual** `|b − A x|/|b|` (computed by
`solve_iterative_bicg_real_residual`) and the iteration count.

`control_solver_bicg_stop -no` (manual Professional 6.376) now also covers
the breakdown exits (continue with a warning + the current solution),
consistent with its max_iter semantics. Note: the record lookup follows the
legacy pattern (active-index check at 0, read at the current control
index) — for a model whose `control_timestep` uses a non-zero control
index, write `control_solver_bicg_stop <index> -no` with the same index
(the `iface_3d_slip` test documents this).

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

## External features / libraries used

None new — plain C++ and the existing `parallel_sys_routine` machinery.
The symmetry check is one serial pass over the element matrices
(≈ one matvec), run once per solve.

## Hardcoded parameters

- `EPS_dAd = 1e-16`, `EPS_TMP = 1e-16` (breakdown thresholds, unchanged).
- `EPS_P = 1e-10` (diagonal guard, unchanged).
- `EPS_SYMMETRY = 1e-8` (symmetry classification threshold; exact symmetric
  physics measures ~1e-16, genuine asymmetries measure 1e-3..1).
- `MAX_PUKNWN_SYMMETRY = 32` (skip the check — fall back to Bi-CG — for
  exotic dof counts; the dense per-element buffer is `(nnol·npuknwn)²`).
- `max_iter = 10·solve_nlocal` (unchanged).
- `check_error = max(bicg_error·|r0|², bicg_error_minimum)` (unchanged
  formula; `bicg_error` default 1e-10, floor 1e-12).

## Pending refactorings / known limitations

- The consolidation coupling block (materi_stress + groundflow_pressure,
  `materi.cc:619` vs `groundfl.cc:116`) has a sign mismatch between the
  off-diagonal blocks (measured `max_rel_asym` up to 1.0) — the coupled
  system is non-symmetric by construction and takes the Bi-CG path.
- `control_solver_bicg_stop` record-index quirk (see above).
- The staggered fixed point (`0.2315×` 2D, `0.143×` road) is the scheme's;
  the monolithic mixed solve (fix C) or a regularized staggered update
  (fix D) is the definitive repair for the user case.
