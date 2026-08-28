# Tochnog solver findings — three paper lines

Date: 2026-08-28
Status: research agenda (findings detected and documented; repair in progress)

Related artifacts in this repository:
- `ProjectDocs/DIAG-SOLVE-MIXTO.md` — full technical diagnosis (paper-grade, 506 lines)
- `ProjectDocs/manual-developer/group_element_selective_reduced_integration.md` — Q4 SRI (opt-in)
- `ProjectDocs/SEGUIMIENTO-CONVERGENCIA.md` — verification log

---

## 1. What we detected

Working on the Tochnog GNU fork (convergence towards Tochnog Professional), in
the `materi_stress_force` post-processing family, we hit a family of wrong
results. The diagnosis (with quantitative evidence) established:

1. **The "mixed u-sigma" system is not mixed at all.** The stress dofs never
   enter the global matrix. Stress advances through a diagonal lumped update:
   the scheme is **staggered (operator-split)**, and the matrix that the linear
   solver iterates is a **velocity-only, symmetric positive definite (SPD)**
   matrix (measured eigenvalues: hex8 1-element λ in [7.47, 111], cond 15;
   zero negative eigenvalues in every case tested).

2. **The Bi-CG stopping criteria declare success without converging.** Three
   of the four exits are dishonest: a breakdown `|dAd| < 1e-16` returns
   `x = 0`; a stagnation test `|Δerror| < 0.1·check_error` returned
   `RC = 0` with error 4.7e-7 against a requested tolerance of 1e-12.

3. **Flat-residual pathology.** For symmetric structures, the first conjugate
   gradient step lands exactly on a flat residual plane: the identity
   `|r|² |Ād|² = 2 (dᵀĀd)²` held to 17 digits because the load excites
   degenerate eigenvalue pairs (λ = 0.1412 at 95% weight). This explains the
   iteration-0 / stagnation / divergence family observed in 3D.

4. **The staggered scheme's fixed point is not the displacement solution.**
   For a 1-element-in-thickness bending beam in 2D, sweeping the outer
   staggered iterations 1..32 converges the section moment to **0.2315×**
   the analytical value, for both the plain quad4 AND the SRI quad4
   (0.23148× vs 0.23163×). The apparent SRI improvement (0.3125×) is a
   transient of the default 2 outer iterations; the full constitutive shear
   re-enters through the `-Bᵀσ` feedback and cancels the SRI at the fixed
   point.

5. **The linear solver is not the bottleneck.** SuperLU and Bi-CG produce
   byte-identical results (same md5). The limitation is the *scheme*, not the
   solver choice.

6. **Real-world validation.** In a 2D plane-strain road base model (subgrade
   soil + concrete slab), integrating stresses in the concrete gave end and
   midspan moments whose absolute values summed to ≈ 0.23·pL²/8, far below
   pL²/8 — physically impossible (for a fixed-fixed beam,
   |M_end| + M_center = pL²/12 + pL²/24 = pL²/8 is exact statics). The
   quad9 element gives the correct answer; quad9 does not lock.

7. **The quad4 shear locking is the classic parasitic shear** (the bilinear
   shape functions cannot represent pure bending u_y ~ -κx²/2, producing a
   spurious γ_xy). The classical remedy — selective reduced integration
   (Hughes) — was implemented as an opt-in keyword
   (`group_element_selective_reduced_integration`): the element stiffness is
   mathematically correct (verified against the exact C matrix: 3 rigid modes,
   no hourglass, exact patch tests) but the end-to-end result remains limited
   by the staggered scheme.

## 2. How we detected it

1. **MSF L2** (2D stress-force integration): quad4 1-in-thickness bending
   measured mom = 0.231×, she = 0.741× (quad9 exact) → shear locking
   documented.
2. **MSF L3** (3D integration): the loaded cantilever fails (A·b ≈ 0,
   dAd = 3.9e-29; hex27 diverges), so 3D validation had to use prescribed
   deformation → solver gotcha documented.
3. **Q4-SRI lote**: element-level verification perfect (exact stiffness,
   rigid modes, no hourglass, patch tests), but end-to-end only 0.31× →
   the element is fine, something else limits it.
4. **Diagnosis lote**: instrumented the assembled matrix. Found: SPD
   velocity-only matrix, staggered scheme, dishonest stopping criteria,
   flat-residual identity, fixed-point sweep (0.2315×), SuperLU ≡ Bi-CG
   (md5 identity).

## 3. Three paper lines

### Line 1 — "The staggered fixed point is not the solution": solver honesty and conditioning in u-sigma FE codes

- **Scope**: the operator-splitting/staggered scheme and the gap between its
  fixed point and the monolithic displacement solution; dishonest iterative
  stopping criteria that declare success without converging; the flat-residual
  pathology for symmetric structures; implications for FE codes that claim a
  "mixed" formulation but are effectively segregated.
- **Evidence**: eigenvalue spectra, A·b breakdowns, the md5 identity,
  fixed-point sweeps, the real-world road-base case.
- **Targets**: Computer Methods in Applied Mechanics and Engineering /
  IJNME (full paper); Advances in Engineering Software or Finite Elements in
  Analysis and Design (shorter empirical paper).
- **Status**: evidence complete; repair (A)+(B) in progress.

### Line 2 — Sensitivity of solver formulations on the same system

- **Scope**: a systematic comparison on one system: Bi-CG vs SuperLU
  (byte-identical — the star result), honest vs dishonest stopping, CG on the
  SPD matrix, MINRES, and the monolithic mixed solve; conditioning across
  element orders (quad4 vs quad9, hex8 vs hex27).
- **Evidence**: the md5 identity, residual histories, conditioning numbers,
  outer-scheme fixed points per configuration.
- **Targets**: empirical/computational conference (ECCOMAS, COMPLAS) or a
  comparative journal paper.
- **Status**: initial measurements exist; needs the repair to complete the
  comparison table.

### Line 3 — Pedagogical: hand calculations catch rubbish FE results

- **Scope**: a case study on verification by hand calculation. The statics
  check |M_end| + |M_center| = pL²/8 for a fixed-fixed beam, applied to a
  seemingly sophisticated 2D plane-strain road base model (soil + concrete)
  whose integrated moments violate it (0.23·pL²/8). Lesson: complex models can
  produce completely rubbish solutions; a one-line hand calculation detects
  it; the physics check comes before trust.
- **Evidence**: the road-base case, the measured chain (locking → staggered
  fixed point → statics violation), quad9 as the "magic fix" and why.
- **Targets**: International Journal of Engineering Education / European
  Journal of Engineering Education.
- **Status**: case study complete; ready to write.

## 4. Repair lines

| Line | Description | Status |
|------|-------------|--------|
| (A) | Honest stopping criteria in the iterative solver (residual-based; no false-success exits) | IN PROGRESS (lote A+B) |
| (B) | Conjugate Gradient for the SPD system (replacing/augmenting Bi-CG) | IN PROGRESS (lote A+B) |
| (C) | Real monolithic mixed u-sigma solve (MINRES / SuperLU with pivoting) | Proposed — definitive fix, large refactor |
| (D) | Regularization of the staggered scheme so its fixed point matches the true solution | Proposed |
| — | Q4 selective reduced integration (opt-in) | DONE — `group_element_selective_reduced_integration` |

## 5. Open questions

- Is the 0.2315× fixed point a property of the staggered scheme itself, or an
  artifact of the dishonest inner solver polluting the outer iterations?
  (Lote A+B will answer.)
- The hex27 divergence (2.16e+13 after 720 iterations) was not reproduced with
  regular cantilevers; the exact failing mesh from the earlier session is not
  in the repository.
- Moment extraction for meshes with 2+ elements in thickness was
  inconclusive (open question).
