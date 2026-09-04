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

---

## 11. The definitive test: Tochnog Professional comparison (2026-08-28)

The Tochnog Professional binary (version 02-08-2026, from the author's
public Drive, "PublicDennis") was obtained and run on the SAME pathological
models that expose the GNU bug. **The bug does not survive in the
Professional.**

### 11.1 Cantilever, quad4, 1 element in thickness (force7q4.dat)

The exact model of the Professional's own `force7.dat` validation test
(which uses quad9) rebuilt with **quad4** elements (L=100, h=10, 2 elements
along x, 1 in thickness, tip load (fx,fy)=(-1.234,-10.)/unit length).
Statics: N = -12.34, V = +100, M(x=50) = -5000, M(x=0) = -10000, M(x=100) = 0.

| quantity | node_dof_calcul value (Professional) | statics |
|---|---|---|
| nory_sig (node 2, x=50) | -12.33999999989 | -12.34 EXACT |
| shey_sig (node 2, x=50) | +100.0000000000 | +100 EXACT |
| momy_sig (node 2, x=50) | -4999.999999991 | -5000 EXACT |
| momy_sig (node 1, x=0)  | -10000.00000002 | -10000 EXACT |
| momy_sig (node 3, x=100) | +3.38e-08 | 0 EXACT |

Targets pass (exit 0). The GNU on the same configuration gives 0.2315×
(§5.4) — the locked section moments.

### 11.2 Fixed-fixed beam, uniform load, quad4 1-in-thickness (ffq4.dat)

The user's road-base check: L=100, h=10, 10 quad4 elements along x, 1 in
thickness, uniform p=1/unit length on the top edge, both ends fixed.
The statics identity |M_end| + |M_center| = pL²/8 = 1250 must hold.

| node | momy_sig (Professional) | expected |
|---|---|---|
| 1 (x=0, end) | -825.0000000036 | -pL²/12 = -833.33 (Euler) / -825 (deep beam, Timoshenko) |
| 6 (x=50, center) | +425.0000000068 | +pL²/24 = +416.67 (Euler) / +425 (deep beam) |
| sum | 1250.0000000104 | pL²/8 = 1250 EXACT |

**|M_end| + |M_center| = 1250 = pL²/8 exactly** (to 1e-10). The individual
values (825/425 instead of Euler 833/417) are the correct deep-beam
(Timoshenko) fixed-end moments for h/L = 0.1 — shear deformation is
physical here, not an error. The statics identity that the GNU violates by
a factor ≈ 7 (§8.1 road model) holds exactly in the Professional.

Nodal shear (shey_sig) shows the same free-surface pollution band seen in
the GNU MSF family (45 vs the exact reaction 50 at the constrained ends;
5 vs 0 at the center) — a shared characteristic of nodal section-force
averaging, unrelated to the lock.

### 11.3 What this proves

- The quad4 1-element-in-thickness locking + staggered-scheme fixed point +
  dishonest solver stopping is **specific to the GNU open-source line**
  (verified back to the 2009/2011-era binary architecture: "stresses follow
  from the principal unknowns", 2011 manual §"initia").
- The Professional produces exact section statics with the SAME element,
  SAME mesh, SAME configuration → its formulation/solver does not suffer the
  bug (consistent with a monolithic mixed solve or an honest solver; the
  exact-to-1e-10 results suggest a formulation whose section equilibrium is
  built in).
- The Professional's own validation tests use quad9 (force7/8) and 2×2
  elements in section for hex8 (force10/13) — a defensive choice that
  sidesteps the GNU-only failure mode; the quad4 1-in-thickness case itself
  works correctly in the Professional.

### 11.4 Reproducibility

- Professional binary: `tochnog_version_02-08-2026` from the author's public
  Drive folder (PublicDennis; tochnog_linux_64_bit.tar.gz, build date
  2-8-2026; statically linked, with debug_info). Not committed to this repo
  (56 MB).
- Test inputs: `/tmp/opencode/tn_prof/force7q4/force7q4.dat` (11.1) and
  `/tmp/opencode/tn_prof/fixedfixed/ffq4.dat` (11.2) — regenerable; the
  Professional syntax is documented in its own `test/other/force*.dat`.
- Professional's own validation family: `test/other/force7.dat` (quad9 2D),
  `force10.dat` (hex8 3D), `force13.dat` (hex8 3D), `beam2d_3.dat`,
  `test/tutorial/tutorial_4/tutorial_4.dat` — all use
  `post_calcul -materi_stress -force` with statics targets.
- **Repeatable baseline**: the full A/B harness (every model in both
  binaries, node_dof_calcul comparison, ratios) and the baseline tables
  live in `ProjectDocs/VALIDACION-PROFESIONAL.md` +
  `scripts/compare_professional.sh` (2026-08-28).


---

## 12. Fix C/D implemented (2026-08-28) — the exact mechanism and the fix

**Status**: DONE (fix D). The staggered scheme's fixed point is now the
element solution; the arness (Professional comparison) is the acceptance
criterion. The monolithic mixed solve (option C) was NOT needed.

### 12.1 The exact mechanism (measured, with code lines)

The staggered loop, iteration by iteration:

```
1. element_loop (elem.cc:404):
   - momentum matrix  dt*K_uu  with  K_uu = B^T*D_elem*B
     (materi.cc:519 matrix_atba + materi.cc:651; SRI: D split, the
     shear integrated at 1 Gauss point at the centroid, materi.cc:654)
   - momentum RHS    P - V*B^T*sigma  with sigma = new_sig, the
     CONSTITUTIVE stress of the current velocity iterate
     (stress.cc:1182-1183: new_sig = sigma_old + C:inc_epe with
     inc_ept = B*v*dt from set_deften_etc, materi.cc:1131;
     materi.cc:508 matrix_atb(new_b, sigvec, force) + materi.cc:617)
2. solve (so.cc:740):  v_new += dv = v + (dt*K_uu)^-1*(P - B^T*sigma)
   (the velocity ACCUMULATES across the equilibrium iterations)
3. parallel_new_dof_diagonal (dof.cc:155): the sigma dofs (never in
   the global matrix: dof_principal unassigned, input.cc:439-444;
   so.cc:298 test2) are updated by the lumped "inertia" equation:
     sigma_new += RHS_sigma/LHS_sigma
   with RHS_sigma = V*h*(new_sig - sigma_old)/dt  (materi.cc)
                   - V*h*(sigma_iter - sigma_old)/dt  (general.cc:253-255,
                   the inertia, ALWAYS active for non-principal dofs)
   and LHS_sigma = V*h/dt (general.cc:256)
   ->  sigma_new = new_sig = sigma_constit(v_iterate)   (REPLACEMENT)
```

The measured dynamics (verified against the dumps and direct solves):

- **The velocity converges in one pass** for the non-SRI elements:
  `v(k+1) = v(k) + (dt*K_uu)^-1*(P - B^T*sigma_constit(v(k)))` with
  `sigma_constit(v) = sigma_old + dt*C_full*B*v` collapses algebraically
  to `v* = (dt*K_uu)^-1*(P - B^T*sigma_old)` when `K_uu = B^T*C_full*B`
  (the plain quad4/quad9/hex8). The dump of the FIRST-solve system
  solved directly equals the run's velocities to the last digit.
- **The fixed point is governed by the FULL constitutive, NOT by the
  momentum matrix.** At the fixed point `dv = 0`, so
  `P = B^T*sigma_constit(v*) = B^T*sigma_old + dt*B^T*C_full*B*v*`
  and the momentum matrix K_uu drops out of the equation:

```
v* = (dt*K_full)^-1 * (P - B^T*sigma_old)        (K_full = B^T*C_full*B)
sigma* = sigma_old + C_full*eps(u*)             (the full-law stress)
```

  The SRI matrix in the momentum equation only shapes the TRANSIENT
  (iteration 2 = the old "0.3125x" SRI value); at equilibrium the
  full-rule shear re-enters through `-B^T*sigma` and the fixed point is
  the LOCKED state `K_full^-1*P` — the SRI benefit is cancelled
  (measured 0.2316x ~ the plain 0.2315x).

- **The sigma dof recovery is lumped-inconsistent for the interior
  quadratures.** The sigma dofs are advanced by the h-weighted average
  `sigma_node = sum_gp V*h*sigma_gp / sum_gp V*h`. With the
  node-containing quadratures (the default 2x2 Lobatto corners of the
  quad4, the quad9/hex8 Lobatto rules) h is the Kronecker delta and the
  recovery is exact. With the 2x2 GAUSS rule (interior points at
  +-1/sqrt(3), switched by the SRI quad4, polynom.cc:421) the h-weighted
  average dilutes the corner values (measured 0.577x of the exact value
  for the bilinear): the section moments read systematically low nodal
  stresses even when the displacement field is correct.

- **The arness divergences (gforce7q4, gforce10/13) are INPUT BC bugs,
  not scheme bugs.** The GNU conversions clamped the WRONG edges:
  - `gforce7q4.dat`: `-ra 1 2 3 -velx` bounded the BOTTOM edge (nodes
    1,2,3 at y=0) instead of the LEFT edge (nodes 1,4) of the
    Professional's `-left_edge`. The missing constraint leaves the
    RIGID ROTATION of the whole cantilever about the clamp node (0,0)
    as a zero-energy mode: the momentum matrix is singular (rank 7/8,
    measured pivot ratio 8.7e-10) and the load projects 88% onto it ->
    the honest CG diverges (residual 1.8e-4 -> 1.4e12). The
    Professional's own test fixes velx on the FULL left edge.
  - `gforce10.dat`/`gforce13.dat`: `-ra 1 4` is the node LIST {1,4}
    (the -ra range_expand processes each integer individually), so the
    bottom face was constrained at only TWO opposite corners and the
    rigid rotation about their diagonal stayed free (singular 30x30
    matrix, pivot ratio 8.6e-17). The correct list is {1,2,3,4} (the
    whole bottom face = the Professional's `-bottom_edge`).

### 12.2 The fix (D: element-consistent feedback + consistent sigma recovery)

**D-c — element-consistent momentum feedback (materi.cc).** The
momentum right-hand side must carry the ELEMENT internal force
`B^T*sigma_old + dt*K_elem*v`, not the full-constitutive stress
`B^T*sigma_constit(v)`. For the SRI quad4 (the only element whose
matrix differs from B^T*C_full*B):

- the current-iterate shear increment of the feedback stress is zeroed
  (`sigvec[stress_indx(0,1)] -= 2*sri_g*inc_ept[1]`, keeping the old
  shear prestress sigma_old_xy), and
- the reduced 1-point shear internal force `-dt*K_shear*v` is added to
  the momentum RHS (the same matrix-vector product the SRI shear block
  adds to the momentum matrix).

With the element-consistent feedback the velocity map collapses to
`v(k) = v* = (dt*K_elem)^-1*(P - B^T*sigma_old)` for every k >= 1:
ONE-pass convergence, no drift, and the fixed point is the element
solution:

```
AFTER:  v* = (dt*K_elem)^-1 * (P - B^T*sigma_old)     (the element)
        sigma* = sigma_old + C_full*eps(u*)           (the physical)
```

For the plain quad4/quad9/hex8 (K_elem = K_full) the feedback is
unchanged: the fixed point and the transients are BYTE-IDENTICAL to the
old scheme (verified on gforce7q4_ref/gffq4/gforce7_ref/qsri_beam2d).

**D-b — consistent sigma recovery (materi.cc + general.cc).** The
sigma dof update uses the bilinear Lagrange EXTRAPOLATION of the
Gauss-point values to the nodes (the "same B at the node") for the
NORMAL stress components, which are superconvergent at the Gauss points
(measured: the SRI gp sigma_xx = 96% of the analytic beam stress; the
extrapolated nodal sigma_xx = 93.9%). The SHEAR components keep the
h-weighting: the Q4 shear is NOT superconvergent (the interpolation
error dominates; the centroid-biased average is the better estimate).
For every node-containing quadrature the extrapolation reduces to h
(the Kronecker delta), so ONLY the SRI quad4 (2x2 Gauss) changes.
The momentum feedback does not read these dofs (it uses the fresh
constitutive stress), so the recovery change is purely an OUTPUT
improvement: the section forces now read the accurate nodal stresses.

### 12.3 Results (arness = acceptance, suite = regression)

- **SRI quad4 cantilever** (qsri_beam2d_sri): the clamp moment goes
  from the 0.2315x locked fixed point to **0.0750 = 0.9375*P*L** (the
  classic Hughes SRI reference; the msf moment about the element
  centroid = P*(L-0.5) exactly), stable at 32 equilibrium iterations
  (no drift). The normal stress field matches the analytic bending
  stress to 94%.
- **Plain quad4**: byte-identical (the locked element solution — the
  lock is ELEMENT physics, opt-in SRI; the scheme now converges to the
  element's own solution).
- **quad9** (gforce7/gforce7_ref): unchanged (M 0.996x/0.9986x).
- **gforce7q4** (2 quad4, BC-fixed): CONVERGES (was diverging); the
  2-element mesh's own locked solution (mom 0.033x — the very coarse
  mesh; documented).
- **gforce10/gforce13** (hex8 3D, BC-fixed): CONVERGE (was diverging);
  the axial N = 12.34 EXACT (1.0000x vs the Professional).
- **Suite**: 199/199 runs + all file checks green; the
  `qsri_beam2d_sri` file check was updated to the fixed 0.9375x moment
  (the old 0.3125x value was the TRANSIENT, which the honest fixed
  point supersedes).

### 12.4 What the fix does NOT change (documented limitations)

- The plain quad4/hex8 shear LOCK (element physics; the opt-in SRI
  keyword fixes the quad4 and, since 2026-08-29, the hex8: the loaded
  8×1×1 cantilever recovers 0.897× of the Euler-Bernoulli deflection
  vs 0.221× locked — with the measured caveat that the shear-only SRI
  hex8 retains zero-energy twist/warping modes, see
  manual-developer/group_element_selective_reduced_integration.md
  HEX8 extension). The scheme now converges to the element's OWN
  solution, locked or not.
- The section SHEAR pollution (the raw sigma_xy of the Q4 — the
  interpolation error; documented in VALIDACION-PROFESIONAL §3 as the
  "mixed shear pollution band"). The Professional's NODAL sigma_xy is
  equally polluted (force7q4 clamp sigma_xy = -90.9 vs the beam
  tau(y=0) = 0); its section statics are exact because they do NOT
  come from the raw nodal stress (consistent with an
  equilibrium/internal-force-based section calculation or a monolithic
  mixed solve with sigma as global unknowns).
- The coarse-mesh lock of the 2-element models (their own formulation
  solution).
- The monolithic mixed solve (option C) is NOT needed: the staggered
  scheme with the element-consistent feedback converges to the element
  solution in one pass.

---

## 13. Fix E (2026-08-29) — the equilibrated sigma state: ELEMENT_DOF populated + section kinematics corrected

**Status**: DONE. Two independent mechanisms were mislabeled as "the
staggered scheme's non-equilibrated sigma state" (VALIDACION §11, the
L4/L5 open front): (a) the ELEMENT_DOF stress block of the 3D
`derivatives` models was never written (the "inc_ept=0" mystery of L4 —
measured mechanism below), and (b) the LOT-5 section internal forces of
the 2D quad9 were computed with a transposed-inverse-Jacobian bug that
amplified the statics by a mesh-dependent factor (5× for the 50×10
gforce7 elements). With both fixed, the recovered sigma satisfies
Bᵀσ = P at the end of the step for every family in the arness: the
"non-equilibrated state" was a post-processing artifact, NOT a scheme
property. The scheme's fixed point WAS the equilibrium all along.

### 13.1 The ELEMENT_DOF zero-stress mechanism (the L4 "inc_ept=0" mystery, resolved)

Measured with per-iteration dumps (element_dof / new_unknowns /
new_dof) on gforce10 (hex8 3D, `derivatives`, nder=5, 1 step × 2
equilibrium iterations):

1. The constitutive stress IS computed at every integration point
   (new_sig = σ_old + C:inc_ept, set_stress) and it IS used — in the
   momentum feedback (materi.cc:559 `matrix_atb(new_b, sigvec, force)`
   and the RHS at materi.cc:717) and in the nodal recovery
   (dof.cc:155, the lumped σ dof update). The run's node_dof σ is
   correct (node 5's σxz = 0.7407 = E·εxz — verified).
2. **The ELEMENT_DOF write never saw it** (materi.cc:850): the stress
   was written to `new_unknowns[stres_indx/nder + j]` — the
   UNKNOWN-NUMBER index — which for nder>1 lands inside the
   displacement block (gforce10: stres_indx/5 = 15 → slots 15..20 =
   disx..disz), NOT at the value slots `stres_indx + j*nder`
   (75..105). The elem.cc write (which copies the stres slots) then
   stored the INTERPOLATED NODAL σ — zero at every write because the
   nodal recovery runs AFTER the element loop in the same iteration
   (the one-iteration lag: the recovery at it=1 reads an RHS built
   with v=0, so the it=2 element loop still sees zero nodal σ, and the
   it=2 recovery — the final one — happens after the last write).
3. A secondary corruption: the restore/write ranges for the
   epe/epp/epi blocks used the uninitialized globals epe_indx=-1 etc.
   when those initia were absent, so `i>=epe_indx && i<(epe_indx+9)`
   matched slots [0,8) and copied the VELOCITY values into the
   ELEMENT_DOF head (the "3.27e-16 0.000377" noise seen in the L4-era
   .dbs). The stress block [stres_indx, stres_indx+6*nder) stayed
   zero, which triggered the L4 fallback ("the element integration
   point stresses are all zero") and its interpretation "the staggered
   loop does not propagate the deformation" — the deformation WAS
   propagated (the RHS and the recovery prove it); only the
   ELEMENT_DOF RECORD was blind.

**The fix (minimal, nder=1 byte-identical)**:
- materi.cc:850: `new_unknowns[stres_indx + (ipuknwn - stres_indx/nder)*nder]`
  (the value slot of the j-th stress component; for nder=1 this is
  exactly the old expression). Same one-line fix for the
  materi_history_variables write (materi.cc:807).
- elem.cc: the stress restore range widened from MDIM*MDIM=9 to
  6*nder (the L4 write already used 6*nder — the restore was
  inconsistent and read a partial/mixed block for nder>1); the hisv
  restore/write ranges now use materi_history_variables*nder (the old
  narrow ranges lost the value slots of the derivatives models) with
  a `hisv_indx>=0` guard.

**Why no extra pass is needed (the one-pass property)**: the scheme
converges in ONE pass (v(k) = v* = (dt·K_elem)^-1·(P − Bᵀσ_old) for
k ≥ 1, fix C/D), so the material call of the LAST element loop runs on
the converged velocity and the element_dof now carries σ* = the
equilibrium constitutive stress. The Bᵀσ = P condition holds for the
written state at the end of the step.

### 13.2 The 5× section statics mechanism (the L5 "non-equilibrated gforce7", resolved)

The arness measured gforce7 (2 quad9, 50×10) at N/V/M = 5.0000× and
gforce7_ref (8 quad9, 12.5×10) at 1.2500×. The σ field is NOT
non-equilibrated: the momentum residual of the last solve is 2.19e-13
(the honest CG, fix A), the velocity is the equilibrium Timoshenko
field (vely(x=50) = -0.01154 vs the analytic 0.01138), and the nodal
σxx(x=50) = ±301 matches the analytic bending stress ±300. The
amplification is a KINEMATICS BUG of the LOT-5 section integration
(msf_element_internal_forces_2d, calcul_force.cc):

- The physical derivatives were computed with the inverse of the
  TRANSPOSED Jacobian: `dn[0] = invjac[0]*pξ + invjac[2]*pη` uses
  invjac[2] (= ∂ξ/∂y) where the correct dN/dx needs invjac[1]
  (= ∂η/∂x); the cross terms invjac[1]↔invjac[2] were swapped.
- The tochnog quad9 local ordering [9,8,7,6,5,4,3,2,1] runs the local
  ξ-axis VERTICALLY (local node 0 = the top-right corner), so the
  Jacobian is the non-diagonal [[0,-5],[-25,0]] and the swap matters:
  dN_0/dx at the corner reads 0.3 instead of the true 0.06 — a 5×
  error for the 50×10 element (the factor = the element aspect ratio
  times the swap; measured 5× for 50/10 and 5/4 for 12.5/10).
- The quad4/quad9-square elements have diagonal Jacobians (the swap
  multiplies zeros) — gforce7q4, gffq4, msf_beam2d (1×1 quad9s) were
  exact and are BYTE-IDENTICAL after the fix. The 3D version uses the
  full matrix-vector product invjac·p3 (correct) — msf_tunnel3d /
  msf_cant3d_hex27 were exact and unchanged.

**The fix (minimal)**: swap invjac[1] and invjac[2] in the two dn
lines of msf_element_internal_forces_2d (the SRI branch's b_shear —
written later — already used the correct formula).

### 13.3 Results (arness = acceptance, suite = regression)

| model | BEFORE | AFTER | Professional |
|---|---|---|---|
| gforce7 (2 quad9) | N/V/M 5.0000× | **N 1.0000×, V 1.0000×, M 0.9984×** | EXACTO |
| gforce7_ref (8 quad9) | N/M/V 1.2500× | **N 1.0000×, V 1.0000×, M 0.9987×** | EXACTO |
| gforce7q4 (2 quad4) | N/V 1.0000×, M 0.9984× | byte-idéntico | EXACTO |
| gffq4 (10 quad4) | identidad pL²/8 0.9996 | byte-idéntico | EXACTO |
| gforce10/13 (hex8 3D) | N 1.0000× (fallback) | **N 1.0000× (ELEMENT_DOF real, sin fallback)** | EXACTO |

The equilibrium check (assembled internal forces from the ELEMENT_DOF
σ with the corrected kinematics, gforce7): Σ f_elem = -P at the free
dofs and the clamp reactions (12.34, 100) match the loads exactly —
Bᵀσ = P to the solver tolerance (1e-5). The L5's "the state σ is not
in equilibrium" (VALIDACION §11) is retracted: the state was always
equilibrated; the L5 kinematics amplified it.

### 13.4 What this does NOT change

- The 3D section V/mom of the axial-loaded arness models (gforce10/13:
  the section shear/moment stay ≈ 0 — a remaining 3D face/arm issue of
  the section post-processing, NOT the scheme; the N is exact and is
  the arness acceptance). The y-loaded 3D cantilever (msf_cant3d_hex27)
  DOES produce the free-body statics (shes = P within 4%).
- The plain quad4/hex8 lock, the SRI hex8 zero modes, the section
  shear pollution of the raw Q4 σ_xy (all documented in §12.4).
- The targets of mesh_act_grav / cmat_gate / qsri3d_beam were UPDATED
  with justification: they were calibrated against the partial σ
  restore of the derivatives models (the old restore read only 9 slots
  and lost σyy/σxy... for nder=4; the fixed restore reads the full
  6*nder block). mesh_act_grav's velx flips from the artifact -0.95 to
  the physical +1.3333 (the +x gravity ramp: the elastic displacement
  u = F·L/(E·A) = 1000·1/1000 = 1); cmat_gate's sigxy from -88.39 to
  +9.09 (the linear shear with the capped first step); qsri3d_beam's
  section moments now read the free-body statics P·L = 0.08 for both
  SRI and OFF (the SRI discriminator is the deflection, not the
  section moment).

---

## 14. Fix F (2026-09-04) — the interface element MEASURE (corpus patch1)

**Status**: DONE. The staggered scheme's fixed point for the interface
elements was the element solution of a system whose assembled interface
forces/stiffness were missing the element LENGTH (2D) / AREA (3D).

### 14.1 The measured mechanism (patch1 of the corpus)

patch1 = two inclined quad6 interfaces (length 1.677 and 0.559) between
two loaded quad9 blocks. The Professional gives sigma_n = 960 uniform per
intpnt (the rotated traction of the sigma_xx = 1200 patch); the GNU
converged (any solver: Bi-CG ≡ band LU ≡ SuperLU, 1-40 equilibrium
iterations, 1 or 2 steps — identical) to sigma_n = (1610, 1610, 1073,
536, 536) per intpnt, mean 1431.08, with the whole right block shifted
~1e-7 (rotated ~7e-8 rad): the σn record = kn·(jump of the elastic
field), and the elastic field itself was wrong.

The root: the assembled nodal force of a closed interface pair was
`w_i * sigma` and the pair stiffness `w_i * kn` with the Lobatto/even
weights w_i (sum = 1) — the integral over the UNIT-length element. The
physical integral is `w_i * L * sigma` / `w_i * L * kn` with L = the
element length (2D) or the face area (3D). The unit-length validations
(interface1-15, conspr1-7, interface_quad4_hex8, ...) never
discriminated the missing measure (L = 1); patch1 did. Without the L the
discrete force system of the interface loses the load-path moment arm:
the length-weighted centroid of the pairs is y = 0 (the load's line of
action) only WITH the measure — the bare Lobatto weights of the unequal
elements 5/6 give a spurious centroid offset, so a uniform traction
cannot balance the applied edge load and the equilibrium sigma_n becomes
the non-uniform pattern above (its weighted sum balances the 2400 load
exactly — verified numerically). With the measure the uniform traction
960 equilibrates (the per-IP sigma_n = 959.99999999-960.00000001, the
node displacements match the Professional's to the 10th digit).

The 3D face area must use a triangle fan (triangle 0-1-2 + triangle
0-2-3): the tochnog hex8/quad4 face numbering (e.g. 5,6,7,8 of the unit
hex8) is a bowtie order — the crossed-diagonals formula gives 0 (the
diagonals coincide). The 2D length = the chord between the first and the
last side-1 node of the reference geometry (the same memory branch as the
interface frame: -total_linear → NODE_START_REFINED).

### 14.2 The solver residual (patch1's ±1e-3 targets)

With the measure the residual of the patch1's sigma_n is the Bi-CG
accuracy on the penalty-conditioned mixed system (kn/E = 1e4): the old
relative criterion (bicg_error 1e-10 → |r|/|r0| < 1e-5) left the
interface jumps noisy at ~1e-2 (sigma_n 960.01). The defaults were
tightened (bicg_error 1e-14, floor 1e-16): the penalty systems now fail
honestly (the CG breakdown — the flat-residual-plane identity of the
symmetric structures) and the existing Bi-CG → direct-LU retry resolves
them exactly (sigma_n = 960.0000000). Well-conditioned systems converge
tighter without the retry; the groundflow large models (ground19) still
converge (18.8 s, no retry). Corpus: 153 PASS (patch1 in; the only
per-test change vs the 141-baseline: validation_2 now fails honestly at
-0.688985 — the same value with the pure direct solver — its -0.63
target was a loose-CG trajectory artifact of the inertia/consolidation
transient; taylor3 = the 45-s timeout marginal, flaky in either binary).

### 14.3 What this does NOT change

- The mpc3/4 tying (0.299/0.350 vs 1/3), ground8 (phreatic_multiple +
  mechanics), dynamic1/2/5/8 (materi_dynamic), ground15/16, mpc5/6:
  different mechanisms (the mpc generation ties the velocity dofs but
  the mixed σ-dofs of the tied nodes stay free → non-homogeneous field;
  the ground8 = the phreatic-multiple/mechanics coupling; the
  dynamic* = the explicit limit of the staggered scheme) — PENDIENTE
  with the fine diagnosis of the next lote.
