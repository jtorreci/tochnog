# DIAG-SOLVE-MIXTO.md — Root-cause diagnosis of the GNU mixed u-σ solve

**Status**: completed (2026-08-28) — branch `documentation-improvement`.
**Fixes A+B implemented** (2026-08-28, §8.1): honest stopping criteria +
CG for the SPD system; the staggered fixed point is confirmed scheme-owned
(0.2315× unchanged) → C/D pending.
**Scope**: the `materi_stress` mixed formulation of the GNU (sfnet 2014
fork): why the default Bi-CG solve produces `A·b ≈ 0` / divergence in 3D,
ill-conditioning and wrong section moments in 2D, and how the user-level
symptom (`|M_end| + |M_center| ≪ pL²/8` in a road-layer model) follows
from the scheme — not from the element shape functions or the integration
of `post_calcul -materi_stress -force`.
**Method**: instrumentation of the assembled system (sparse/dense matrix
dumps + RHS + dof maps, gated by the `TOCHNOG_DUMP` env var, removed
afterwards), analysis with numpy (spectra, `A·b`, quadratic forms,
null-space projections), A/B solver comparisons (Bi-CG vs SuperLU) and
equilibrium-iteration sweeps. Every number below was measured on the
current build; the diagnostic scratch files live in `/tmp/opencode/diag`.

---

## 1. Executive summary

The GNU "mixed u-σ" solve is **not** a monolithic saddle-point system:
with `materi_stress`, the stress dofs are *never part of the global
matrix* (`dof_principal` is never assigned to them in `input.cc`), so the
Bi-CG system is a **velocity-only, symmetric positive definite (SPD)
system** (verified: all eigenvalues positive in every dumped case). The
stress field is updated *after* each solve by a lumped-diagonal
"constitutive update" — the scheme is **staggered**, not mixed.

Two independent, measured failure mechanisms follow from this design:

1. **The Bi-CG stopping criteria are not convergence criteria.** Three of
   the four loop exits report *success* without requiring
   `error < check_error`:
   - `dAd < 1e-16` (or `r1·r2 < 1e-16`) → **stop at iteration 0 with
     x = 0** whenever the load vector is (numerically) in a
     (near-)null mode of the preconditioned matrix;
   - `|error − last_error| < 0.1·check_error` (stagnation) → **stop with a
     partial solution whose L2 residual is far above the tolerance**.
   The L2 residual of CG/Bi-CG is *non-monotone*, so a flat or rising
   step is normal (measured: errors rising over the first iterations of
   well-posed 3D cases), and for mirror-symmetric structures the first CG
   step can land *exactly* on a flat residual (measured identity
   `|r|²·|Ād|² = 2·(dᵀĀd)²` to 17 digits), which triggers the
   stagnation exit with a wrong solution reported as `RC=0`.
2. **The staggered scheme's fixed point is not the displacement
   solution.** The momentum matrix contains the full stiffness
   `dt·K = dt·Bᵀ·C·B` *and* the right-hand side contains the stress
   gradient `−Bᵀ·σ`, with σ = C:ε(v)·dt at equilibrium — the stiffness
   operator is effectively applied twice, and the constitutive shear
   re-enters through σ even when the element matrix was selectively
   reduced (SRI). Measured on the 1-element-in-thickness quad4
   cantilever: the fixed point is `mom = 0.2315×` **for both the plain
   and the SRI element** (SRI benefit cancelled), the documented
   `0.3125×` SRI value is a transient of the default 2 equilibrium
   iterations, and **SuperLU produces byte-identical results to Bi-CG**
   (md5 `62e82a92...`) — the linear solver is not the limiter, the
   scheme is.

**User case**: a 1-element-in-thickness quad4/hex8 layer under bending
locks, and the staggered scheme converges to the locked state
(`≈ 0.23×` the section moments). For a fixed-fixed slab the statics
requires `|M_end| + |M_center| = pL²/8` exactly (derivation in §7); the
FE result `≈ 0.23·pL²/8` violates it by a factor ≈ 4 — the static check
the user applied is exactly the right detector.

---

## 2. Symptom

### 2.1 Documented in the MSF sub-sprint (GOTCHA of lot 3)

- **3D hex8 / hex27 cantilever, real tip load**: the default solve
  "gives `A·b ≈ 0` and the Bi-CG stops at iteration 0 with x = 0"
  (measured at the time: `dAd = 3.9e-29`), or diverges (hex27:
  `final error 2.16e+13` after 720 iterations = `10 × solve_nlocal`).
- **3D ring under pressure**: converges poorly (plateau `~1e-12`),
  contaminated stresses.
- **Consequence**: the 3D MSF tests use *prescribed deformation*
  (Dirichlet, no solve); load-based validation is blocked.

### 2.2 Documented in the 2D MSF/SRI work

- 1-element-in-thickness quad4: the mixed system is "ill-conditioned"
  (`cond ≈ 2·10⁴`); Bi-CG gives moments `0.31×` the analytic value,
  while "exact Gaussian elimination" was reported to give `93.75%`
  (claimed: the SRI element is correct, the solver is the limiter).
- quad9 (no lock): the same Bi-CG converges perfectly.

### 2.3 User case (road layers, 2D plane strain)

Soil + concrete layers; section forces integrated over the concrete
layer via `post_calcul -materi_stress -force`. Result:
`|M_end| + |M_center| ≪ pL²/8` — statically impossible. The user
attributed it to the element formulation; the shape functions and the
integration are fine — the interaction element↔solver is the problem.

---

## 3. Minimal reproduction

All cases below: `E = 1000, ν = 0.3`, cantilever `L=8, h=1 (b=1 in
3D)`, clamped at `x=0`, total tip load `P = 1e-2` in −y
(`bounda_force -vely`, consistent face distribution), 1 equilibrium
iteration (`control_timestep_iterations 1 1`) unless stated.

### 3.1 3D — 1 hex8 element, tip load (stagnation with wrong solution)

```
echo -no
number_of_space_dimensions 3
derivatives
materi_velocity
materi_velocity_integrated
materi_strain_total
materi_stress
number_of_integration_points 8
end_initia
options_mesh -fixed_in_space -fixed_in_space -fixed_in_space
options_convection -no
options_inertia -no
node 1 0 0 0   ... (unit cube, 8 nodes)
element 1 -hex8 1 2 3 4 5 6 7 8
group_type 0 -materi
group_materi_elasti_young 0 1000.0
group_materi_elasti_poisson 0 0.3
bounda_unknown 0 1 -velx -vely -velz   (4 clamped nodes, x=0)
bounda_force 4 2 -vely   (4 tip nodes, 2.5e-3 each)
control_timestep_iterations 1 1
control_timestep 1 0.1 0.1
end_data
```

Result: `RC=0`, solver prints `initial error 4.72727e-07`,
`final error 4.72727e-07` (check error is `1e-12`) — **the run reports
success without solving**. With `print_solver -yes`:

```
0 4.72727e-07
1 4.72727e-07        <- EXACTLY flat; stagnation exit fires
initial error: 4.72727e-07   check error: 1e-12   final error: 4.72727e-07
```

### 3.2 2D — 8×1 quad4 with SRI, distributed tip load

The quad4-SRI matrix is **rank-sufficient** with a correct mesh
(`cond 2.04e4`, full rank); with an **invalid (counter-clockwise) node
order** the matrix loses rank (`rank 24/32`, 8 null modes) and the load
falls exactly in the null space: `A·b ≈ 1.4e-14·|b|`, `dAd = 0`, Bi-CG
stops at iteration 0 with `x = 0` (`dAd = -8.7e-23` preconditioned) —
the exact `A·b ≈ 0` breakdown family. Lesson: tochnog quad4/hex8
connectivity must follow the **Z convention** (`1 2 10 11`, bottom row
then top row in x-order); CCW ordering silently produces negative
Jacobians and singular matrices.

### 3.3 What works

| Case | Bi-CG iterations | final error |
|---|---|---|
| 3D hex8 ×1, top-face load | 6 | 2.07e-36 |
| 3D hex27 ×8, tip load | 79 (error rises first: 4.8e-7 → 1.9e-6 → …) | 6.5e-13 |
| 3D hex8 ×8, tip load | 19 | 1.3e-28 |
| 3D sheet pile 4×hex27, tip load | 22 | 1.9e-13 |
| 2D quad4 ×8 (no SRI) | 12 | 8.4e-13 |
| 2D quad9 ×8 (msf_beam2d) | — | converges, mom 1.03× |

---

## 4. Theoretical background

### 4.1 Saddle-point (mixed) systems

A genuine displacement-stress mixed method leads to the saddle-point
system

```
[ 0   Bᵀ ] [u]   [f]
[ B  −C⁻¹] [σ] = [0]
```

whose matrix is **indefinite** (eigenvalues of both signs); CG-type
methods cannot be applied directly and the standard solvers are MINRES,
GMRES or a direct method with pivoting. Standard reference:
G. Benzi, G. H. Golub, J. Liesen, *"Numerical solution of saddle point
problems"*, Acta Numerica 14 (2005) 1–137.

**The GNU does not assemble this system.** With `materi_stress` the
stress dofs are registered (`input.cc:441`, 6 slots in 2D and 3D) but
`dof_principal` is never set for them (unlike velocity/temperature/
pressure dofs), so `solve()` (`so.cc:298`, `test2 =
dof_principal[iuknwn]>=0`) excludes them from the global system. The
stress is advanced by the lumped-diagonal update
(`parallel_new_dof_diagonal`, `top.cc:414`): `σ_new += RHS_σ/LHS_σ`
with `RHS_σ = V·h·(σ_constit(v) − σ_old)/dt`, `LHS_σ = V·h/dt` — i.e.
`σ_new = σ_constit(v)`, the constitutive stress of the current velocity
field. The momentum equation keeps the full stiffness `dt·Bᵀ·C·B` in the
matrix *and* the stress gradient `−Bᵀ·σ` in the right-hand side
(`materi.cc:609`). The scheme is **staggered**:

```
v_{k+1} = (dt·K)⁻¹·(P − Bᵀ·σ_k),   σ_{k+1} = C:ε(v_{k+1})
```

### 4.2 Bi-CG / CG on an SPD system

For a symmetric matrix the Bi-CG iterates coincide with CG. Two
properties matter here:

- **The L2 residual is non-monotone.** CG minimises the *A-norm* of the
  error; `‖r_k‖₂` can rise (measured: 3D hex27, `4.8e-7 → 1.9e-6` over
  the first iterations). Any convergence test based on a *single-step*
  residual decrease is unsafe.
- **Breakdowns.** The scalar `dAd = dᵀ·A·d` is the denominator of the
  CG step. For an indefinite A it can vanish (true breakdown); for a
  singular or near-singular SPD A with `b` in the (near-)null space it
  is ≈ 0 and the first step is meaningless. The code
  (`so_bicg.cc:120`) treats `|dAd| < 1e-16` as *success* with `x = 0`.

### 4.3 The exact flat-residual identity (new measurement)

For the 1-element hex8 cantilever (and any symmetric structure whose
load excites a degenerate eigenpair — mirror symmetry produces repeated
eigenvalues), the preconditioned quantities satisfy, to machine
precision:

```
|r|² · |Ā·d|² = 2·(dᵀ·Ā·d)²        (measured ratio 1.00000000000000000)
```

with `d = r = P·b` (preconditioned residual) and `Ā = P·A·P`. The CG
step `α = |r|²/(dᵀĀd)` then overshoots the L2-optimal step
`α* = (dᵀĀd)/|Ād|²` by exactly a factor 2, and the quadratic
`‖r − α·Ād‖²` (a parabola in α) returns to its starting value:

```
‖r − α·Ād‖² = ‖r‖² − 2α·dᵀĀd + α²·|Ād|² = ‖r‖² − 2·|r|² + 4·|r|²/2·(...)  = ‖r‖²
```

measured `err1/err0 = 1.0000000000000004` in exact arithmetic. The load
of the symmetric cantilever projects almost entirely onto the degenerate
pair `λ = 0.1412` (95 % of `‖b‖²`) — this is *systematic* for symmetric
structures, not a numerical accident. The stagnation exit then reports
success at `error = 4.7e-7 ≫ 1e-12`.

### 4.4 SRI (Hughes)

Selectively reduced integration for the bilinear quad4: the shear term
`γ_xy` is integrated with one point at the centroid, the normal terms
with the full 2×2 rule. The *isolated* SRI element is correct
(documentation of `group_element_selective_reduced_integration`): 3
rigid modes at zero energy, the 4th (hourglass) mode with positive
energy, and the pure-displacement system `K_SRI·u = P` reproduces the
classic SRI result (tip deflection 99.3 %, clamp moment 93.75 % —
verified by Gaussian elimination on the extracted element matrices).
Reference: T. J. R. Hughes, *The Finite Element Method* (SRI, hourglass
control).

---

## 5. Investigation trail (evidence)

### 5.1 The assembled system is velocity-only and SPD

Dump of the exact matrix the Bi-CG iterates on (element matrices +
`node_lhside` diagonal, rebuilt in `so_bicg.cc` during diagnosis; the
assembly in `so.cc` stores the same data):

| Case | n | eigenvalues (min/max) | negative | near-zero | rank | cond |
|---|---|---|---|---|---|---|
| 3D hex8 ×1, tip load | 12 | 7.47 / 111.1 | 0 | 0 | 12 | 1.5e1 |
| 3D hex27 ×8, tip load | 432 | 1.24e-3 / 7.82e2 | 0 | 0 | 432 | 6.3e5 |
| 2D quad4 ×8 (plain) | 32 | 4.76e-2 / 3.47e2 | 0 | 0 | 32 | 7.3e3 |
| 2D quad4 ×8 (SRI) | 32 | 1.30e-2 / 2.65e2 | 0 | 0 | 32 | 2.0e4 |
| 2D quad4 ×8 (SRI, CCW order) | 32 | ~0 / 4.51e2 | 4 (1e-14) | 8 | 24 | 4.4e17 |

**No negative eigenvalues in any valid case → not a saddle point, not
indefinite.** The "cond ≈ 2·10⁴" of the documentation is the *velocity*
matrix conditioning, which is unremarkable for this mesh family.

### 5.2 RHS and stress-index checks (no assembly bug)

- First-solve RHS dump (1 equilibrium iteration): `b = +0.0025·ŷ` at
  each of the 4 free tip nodes — **the applied load enters the velocity
  dofs only, exactly as specified** (`bounda.cc:779`,
  `node_rhside[ipuknwn] = factor·load`). No leakage to stress dofs.
- `stress_indx(i,j)` (`miscel.cc:321`) is a *dense* 6-component Voigt
  map `(0,0)→0, (0,1)→1, (0,2)→2, (1,1)→3, (1,2)→4, (2,2)→5` — no
  "9-slot with holes"; the momentum gradient (`materi.cc:609`), the
  constitutive RHS (`materi.cc:689`) and the post-processor
  (`calcul_force.cc:1047`) all use the same map. The 2D/3D layouts are
  consistent. (The "holes" hypothesis is **rejected**.)
- Null space of the σ block: **not applicable** — σ is never in the
  global solve; the σ-σ "block" is the positive lumped diagonal
  `V·h/dt` used only by the update.

### 5.3 The breakdowns, measured

| Symptom | Measured evidence | Mechanism |
|---|---|---|
| "Bi-CG stops at iteration 0 with x=0" | 2D SRI + CCW mesh: `dAd = −8.7e-23`, `A·b/‖b‖ = 1.4e-14`, `b` in null space (projection 1.0000000000000002) | `dAd < 1e-16` exit |
| "Bi-CG stops with a wrong solution" | 3D hex8 ×1 tip: error flat at `4.72727e-07` (exact identity, §4.3); `RC=0` with `error ≫ 1e-12` | stagnation exit `|Δerror| < 0.1·check_error` |
| "Slow / unreliable convergence" | 3D hex27 ×8: 79 iterations, error *rising* `4.8e-7 → 1.9e-6` before decreasing | non-monotone L2 residual vs error-based criteria |
| "Divergence 2.16e+13 after 720 iterations" | documented (scratch mesh of the MSF session; not reproduced with the regular cantilevers — same criterion family) | error-based exits at `max_iter` |

### 5.4 The 2D fixed point is not the displacement solution

Section moment at `x = 0` (analytic `P·L = 0.08`), 8×1 quad4
cantilever, SRI on/off, as a function of the equilibrium iterations:

| iterations | plain quad4 | SRI quad4 |
|---|---|---|
| 1 | 0 (σ not yet updated) | 0 |
| 2 | 0.018518 (= 0.23148×, **fixed point already**) | 0.0249998 (= 0.3125×, transient) |
| 4 | 0.018518 | 0.019312 |
| 8 | 0.018518 | 0.018530 |
| 16 / 32 | 0.018518 | 0.018530 (= 0.23163×, fixed point) |

- **The SRI fixed point (0.23163×) equals the plain fixed point
  (0.23148×) within 0.06 %** — the SRI's stiffness correction is
  cancelled at equilibrium: the momentum RHS carries the *full*
  constitutive stress (shear included) through `−Bᵀ·σ`, restoring the
  locked shear stiffness the SRI removed from the matrix.
- The documented "SRI improves 0.231× → 0.312×" is a **transient of the
  default 2-iteration setting**, not the converged state.
- The "93.75 % with Gaussian elimination" is the **pure-displacement
  reference** `K_SRI·u = P` (verified on extracted element matrices in
  the SRI documentation) — *not* the system the GNU assembles.
- **Bi-CG vs SuperLU on the same assembled system: byte-identical**
  section forces (md5 `62e82a927f03c8e933f9a069d842f44a`,
  `mom = 0.0249998465889` at 2 iterations). The linear solver is not the
  differentiator.

### 5.5 Assembly facts (trace of the mixed path)

- `input.cc:441` — `materi_stress` → `stres_indx`, `n = 6` slots,
  `dof_scal_vec_mat = -MATRIX`; **no `dof_principal` assignment**
  (contrast `materi_velocity`, `temp`, `pres`, …).
- `so.cc:293-305` — the global system contains only
  `!node_bounded && dof_principal ≥ 0` dofs → velocity only.
- `materi.cc:609` — momentum RHS stress gradient `−V·Bᵀ·σ`
  (Green's partial integration); `materi.cc:689` — constitutive RHS
  `V·h·(σ_constit − σ_old)/dt`; `general.cc:210` — the σ-σ diagonal
  `V·h·1/dt` (inertia, lumped).
- `top.cc:414` — `parallel_new_dof_diagonal` performs the lumped σ
  update after each solve → staggered loop.
- `so_bicg.cc:102-176` — the four loop exits and their consequences
  (§5.3). `check_error = max(1e-10·‖r₀‖², 1e-12)`.

---

## 6. Root cause (with evidence)

1. **The solve architecture**: `materi_stress` does not produce a
   saddle-point system — the Bi-CG system is the velocity-only SPD
   system and σ is advanced by a lumped staggered update. The
   hypotheses "indefinite spectrum → BiCG unsuitable", "3D assembly
   bug (RHS/indexing)" and "null space of the σ block" are **rejected
   by the dumps** (§5.1–5.2): SPD spectra, correct load RHS, dense
   consistent stress indices, σ absent from the solve.
2. **Bi-CG exit criteria** (`so_bicg.cc`) treat breakdown and
   stagnation as *success*: `|dAd| < 1e-16` → x=0; flat L2 residual →
   "converged" at `error ≫ check_error`. The L2 residual is
   non-monotone, and for symmetric structures the first CG step can
   land exactly on a flat residual (§4.3) — so wrong solutions are
   reported as `RC=0`. This explains the 3D failures (stagnation /
   iteration-0 stop / slow rising-then-falling errors) and the 2D
   "misleading residual".
3. **The staggered fixed point** (§5.4) is not the displacement
   solution: the momentum matrix and the σ-gradient coupling apply the
   stiffness twice and re-introduce the full constitutive shear,
   cancelling SRI. The 1-in-thickness quad4/hex8 then converges to the
   *locked* state (`≈ 0.2315×` moments) regardless of element
   stiffness or linear solver (Bi-CG ≡ SuperLU, byte-identical).

The documented "cond ≈ 2·10⁴ → solver is the limiter" is a
misattribution: the conditioning is unremarkable, the direct solver
gives the same result, and the *scheme* (not the solver) is the
limiter.

---

## 7. Implications

### 7.1 The user's road-layer case (the full chain)

1. **Geometry/formulation**: a 1-element-in-thickness quad4/hex8 layer
   under bending cannot represent the flexural field
   (`u_x ~ y²`); parasitic shear locks the element.
2. **Scheme**: the staggered mixed iteration converges to the locked
   state (§5.4): `mom ≈ 0.23×` the correct value; the SRI keyword does
   not rescue it (cancelled at the fixed point).
3. **Statics check**: for a fixed-fixed beam under uniform load `p`,
   equilibrium gives `M(x) = M_e + (pL/2)x − px²/2`, so
   `M_c − M_e = pL²/8`; with `M_e = −pL²/12`, `M_c = pL²/24` and

   ```
   |M_end| + M_center = pL²/12 + pL²/24 = pL²/8   (exact, fixed-fixed)
   ```

   Hence `|M_e| + |M_c| ≪ pL²/8` (measured `≈ 0.23·pL²/8`) violates the
   moment equilibrium of the section — the element shape functions and
   the MSF integration are exonerated; the element↔solver interaction
   is the culprit, exactly as the user suspected of the wrong component.

### 7.2 Element-level vs end-to-end (SRI)

The SRI element is correct (pure-displacement reference 93.75 %,
documented with Gaussian elimination on extracted element matrices);
the end-to-end result is limited by the *scheme*, not by the element
and not by the linear solver (§5.4: Bi-CG ≡ SuperLU byte-identical).

### 7.3 MSF 3D validation

Load-based 3D validation (`post_calcul -materi_stress -force` with real
loads) remains blocked until the solve is fixed: with the current
criteria a run can report `RC=0` with an unsolved (or garbage) field.
The prescribed-deformation tests of MSF lot 3 remain valid (no solve).

---

## 8. Proposed fixes (with tradeoffs)

**A. Fix the Bi-CG exit criteria (smallest change, global effect).**
Only report success when `error < check_error`; turn stagnation into
"keep iterating (the L2 residual of CG is non-monotone)" and let
`max_iter` fail loudly. Optionally replace the L2 monitor with the
relative residual `‖b − A·x‖/‖b‖` of the *unpreconditioned* system.
Tradeoffs: honest failures instead of silent garbage (tests that today
"pass" via the stagnation path will fail loudly — the 199-run suite
must be re-validated; the flat-residual cases will hit `max_iter` and
fail, which is the desired behaviour for validation but a behaviour
change for users).

**B. Proper solver for the velocity system.** The velocity matrix is
SPD — plain CG (Lanczos) is the natural method; Bi-CG on symmetric
systems reduces to it anyway, so the gain is mainly in clean criteria
(CG's A-norm monotonicity) — it does *not* by itself cure the flat
step. Pair with A.

**C. Monolithic mixed solve (architectural).** Assemble the true
saddle-point system `[0 Bᵀ; B −C⁻¹]` (remove `dt·K` from the momentum
block when `materi_stress` is active) and solve with MINRES (Benzi–
Golub–Liesen) or SuperLU with pivoting (already linked; `-matrix_superlu`
works today). Tradeoffs: correct physics and a well-posed problem; large
refactor of `elem.cc`/`so.cc`; the indefinite system needs MINRES/GMRES
or a direct solver with pivoting; new validation needed for the whole
suite.

**D. Regularisation of the staggered loop.** Add a relaxation/under-
relaxation to the σ update (`σ ← σ_old + ω·(σ_constit − σ_old)`) or a
penalty on the σ-σ block so the fixed point approaches the displacement
solution. Tradeoffs: small change, but the fixed point changes with `ω`
and the equivalence with the true mixed solution is only approximate;
needs careful tuning and re-validation.

**E. Pragmatic, today:** for validation runs use
`control_options_solver -matrix_superlu` (robust pivot handling, no
Bi-CG breakdown) and **more than one element in the thickness** for
bending-dominated sections (removes the lock that feeds the wrong fixed
point). Note that E does not cure the scheme (§5.4: SuperLU = Bi-CG on
the current assembly); it protects against the breakdown family.

Recommended order: A (correctness of the stop criteria) → B (CG) → C
(monolithic mixed, the real fix for the user case) as a separate lot;
D/E as short-term mitigations.

---

## 8.1 Results of lot A+B (implemented 2026-08-28)

**Status**: A (honest criteria) and B (CG) are **DONE** (commit
`feat(solver): ...`); C/D remain the definitive repair for the user case.

### What was implemented

- **A — honest stopping criteria** in `so_bicg.cc`: success is reported
  ONLY on the residual test `error < check_error` (relative, with the
  absolute floor; residual recomputed from the iterate every pass). The
  breakdown exits (`|dAd|<1e-16`, `|r1·r2|<1e-16`) and the stagnation exit
  (`|Δerror|<0.1·check_error`) no longer report success: a genuine
  breakdown is an honest failure (RC≠0, message with the real relative
  residual), a flat step just keeps iterating, and `max_iter` fails loudly
  (with `control_solver_bicg_stop -no` as the documented "continue"
  escape).
- **B — CG for the SPD system**: a runtime symmetry check
  (`solve_iterative_bicg_symmetric`, one serial pass over the element
  matrices, `EPS_SYMMETRY=1e-8`) dispatches to plain CG for the
  measured-symmetric systems (velocity/temp/pres) and to the honest Bi-CG
  for the measured non-symmetric ones (beam: translation↔rotation blocks
  carry different `dtime` factors; plastic-slip 3D interface: `max_rel_asym`
  up to 1.0, `dᵀAd<0`). Same dispatch, same RC contract.
- **Primal-residual monitor**: the accumulated residual is now the PRIMAL
  `b̂−Āx̂` (the old accumulation measured the transpose residual, which
  never vanishes at the solution of a non-symmetric system — the beam
  solve could only "finish" through the false breakdown exit). Identical
  values for symmetric systems.

### The KEY question — did the fixed point change?

**No.** Sweeping the outer staggered iterations 1..32 with the honest
solver gives the SAME fixed point as §5.4:

| iterations | plain quad4 (old) | plain (new) | SRI (old) | SRI (new) |
|---|---|---|---|---|
| 2 | 0.018518 | 0.0185184 | 0.0249998 | 0.0249999 |
| 8 | 0.018518 | 0.0185184 | 0.018530 | 0.0185303 |
| 32 | 0.018518 (0.23148×) | 0.0185184 | 0.018530 (0.23163×) | 0.0185298 (0.23162×) |

The `0.2315×` fixed point is **a property of the staggered scheme**, not
an artifact of the dishonest inner solver (the 2D inner solves were already
converging honestly: `error 8.4e-13 < check_error`). **C/D are still
required.**

### The 3D family is fixed

- 1-hex8 cantilever with tip load (the `A·b≈0` flat-residual case): the
  old code stopped at iteration 1 with `RC=0` and error `4.72727e-07`;
  the new CG iterates through the flat step and converges
  (`iterations 0..3`, `final error 4.2e-37`, `relative residual 9.4e-16`).
- hex8 ×8 and hex27 ×8 with REAL load + MSF section forces: both solve
  (RC=0, CG path). hex27 gives `mom ≈ 1.1·P·(L−x)` and `she ≈ P` — the
  load-based 3D validation of the MSF family is unblocked. hex8 ×8 gives
  `mom ≈ 0.26·P·L` (the 1-in-thickness lock, same family as 2D).
- The road model (concrete 1-in-thickness on soil, plane strain,
  fixed-fixed, distributed load): now runs RC=0 with converged solves, but
  `|M_end| + M_center = 0.001435 = 0.143·pL²/8` — the statics check still
  fails by a factor ≈ 7 (the scheme's fixed point for that model).

### Suite impact (199 tests, clean build)

`198/199` unchanged; `iface_3d_slip` updated (test corregido): the
plastic-slip 3D interface system is non-symmetric/indefinite and the
honest Bi-CG breaks down; the test now opts into the documented
`control_solver_bicg_stop -no` (continue with warning) and keeps
verifying the interface physics. All file checks pass. No test approaches
the 120 s timeout (slowest ≈ 1.4 s). Details: `ProjectDocs/manual-developer/
solve_iterative_bicg.md`.

---

## 9. What each paper line can take from this

1. **Paper line 1 — conditioning of the u-σ solve + solver**: measured
   spectra of the velocity subsystem (SPD; `cond 15 … 6.3e5` across the
   case matrix), the exact flat-residual identity for degenerate
   eigenpairs (§4.3), the failure modes of error-based CG stopping
   criteria, Bi-CG vs SuperLU parity (§5.4). This corrects the
   "indefinite saddle-point" framing: the GNU system is *staggered and
   SPD*, and the interesting mathematics is in the *criteria* and the
   *staggered fixed point*, not in indefiniteness.
2. **Paper line 2 — sensitivity of solver formulations**: the A/B
   table (Bi-CG vs SuperLU byte-identical; iteration sweeps 1→32;
   mesh-order sensitivity with the Z-convention lesson; non-monotone
   residual histories) is a ready-made sensitivity study of the
   solve path for the mixed formulation.
3. **Paper line 3 — pedagogical hand calculation**: the statics check
   `|M_end| + |M_center| = pL²/8` for the fixed-fixed beam (§7.1) is a
   clean, exact, hand-checkable detector that exposes an FE section-
   force error of any origin — an ideal worked example for a
   validation/quality paper.

---

## 10. Files and evidence

- Diagnostic scratch (matrix dumps, spectra, iteration logs):
  `/tmp/opencode/diag/` (not committed; regenerable with the
  `TOCHNOG_DUMP` instrumentation described in §5.1).
- Reproduction inputs: `/tmp/opencode/diag/cant_hex8_1e_tip.dat`,
  `fp2_{plain,sri}_it*.dat`, `ab_qsri_{bicg,superlu}.dat`.
- Code touched during diagnosis: `so.cc`, `so_bicg.cc` (temporary dump,
  **reverted**; `git checkout` — working tree clean, suite 199/199
  green after rebuild).
- Related documentation: `ProjectDocs/manual-developer/
  post_calcul_materi_stress_force.md` (the GOTCHA this document
  explains and corrects), `ProjectDocs/manual-developer/
  group_element_selective_reduced_integration.md` (SRI; the "93.75 %"
  reference is the pure-displacement system), `ProjectDocs/
  SEGUIMIENTO-CONVERGENCIA.md` (tracking row added).
