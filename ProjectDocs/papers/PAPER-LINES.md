# Tochnog solver findings — three paper lines

Date: 2026-08-28
Status: research agenda (findings detected and documented; fixes A+B **DONE**
2026-08-28 — honest stopping criteria + CG; **fix D DONE** 2026-08-28 —
the staggered fixed point is now the element solution, verified against
the Tochnog Professional binary; C not needed)

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
- **Status**: evidence complete; repair (A)+(B) **DONE** (2026-08-28): the
  honest solver converges the 3D flat-residual cases and confirms the
  0.2315× fixed point is the scheme's — the paper's central contrast is
  now "dishonest criteria vs honest criteria vs scheme fixed point" with
  measured numbers for all three. **Repair (E) DONE** (2026-08-29): the
  end-of-step σ state satisfies Bᵀσ = P (the ELEMENT_DOF 3D write fixed
  + the section kinematics corrected) — the last arness gap (gforce7
  5×) closes to 1.0000×, measured as a post-processing artifact, not a
  scheme defect. NEW DATA (sprint 12 lot 1, 2026-08-29): the FULL 3D
  statics of gforce10/13 (nory/shey/mom1y = −12.34/−100/−5000) now
  match the Professional EXACTLY — the equilibrium-based section
  resultants of a converged-but-LOCKED hex8 solution still satisfy
  the exact free-body statics (the lock lives in the stress field and
  the deflection, never in the equilibrium resultants). Paper-1
  implication, measured: "the scheme's fixed point is not the
  displacement solution, BUT its equilibrium resultants are still
  exact statics — what is wrong is the field, and any validation that
  reads statics (instead of deflections) can be fooled". Also: the
  reference Professional binary itself needs `thickness_switch -yes`
  on the square-section cantilever — without it, ITS section frame
  degenerates too (vectors ≈ 0, moment in the wrong slot): validation
  keywords, not physics, decide what a reference implementation
  answers on degenerate input.

### Line 2 — Sensitivity of solver formulations on the same system

- **Scope**: a systematic comparison on one system: Bi-CG vs SuperLU
  (byte-identical — the star result), honest vs dishonest stopping, CG on the
  SPD matrix, MINRES, and the monolithic mixed solve; conditioning across
  element orders (quad4 vs quad9, hex8 vs hex27).
- **Evidence**: the md5 identity, residual histories, conditioning numbers,
  outer-scheme fixed points per configuration. New (2026-08-29): the
  element-order sensitivity row now includes the SRI hex8 — the
  conditioning/zero-mode comparison of the shear-only selective
  integration of the trilinear brick (9 zero modes isolated vs 6 for
  the full rule; the mesh section-warping modes) alongside the measured
  cantilever recovery (0.897× of EB vs 0.221× locked).
- **Targets**: empirical/computational conference (ECCOMAS, COMPLAS) or a
  comparative journal paper.
- **Status**: initial measurements exist; repair (A)+(B) **DONE** adds the
  honest-vs-dishonest and CG-vs-Bi-CG comparison rows (measured); the
  monolithic mixed comparison (C) remains future work. New (2026-08-29):
  the element-order sensitivity row now includes the section-kinematics
  comparison — the LOT-5 transposed-inverse bug amplified the statics
  by the element aspect ratio (5× for 50×10, 5/4× for 12.5×10 quad9s;
  zero for the diagonal Jacobians and the 3D full matrix-vector
  product) — a clean, measured example of post-processing kinematics
  sensitivity on the SAME converged state. NEW DATA (sprint 12 lot 1,
  2026-08-29): a complete measured TAXONOMY of section-frame
  sensitivities, same converged gforce10 state, four independent
  effects: (a) lexicographic border tables → the corner "edges" were
  edge+DIAGONAL (with -yes the diagonal wins the extent rule: a
  45°-rotated frame, nors 17.45 = 12.34·√2 — the error carries √2);
  (b) the square-section tie (10×10) defeats the extent rule — the
  tie-break (reference-point direction vs edge order) decides which
  axis is "thickness", i.e. WHICH slot carries the moment (mom1 5000
  vs mom2 5000, she 100 vs she ≈ 0); (c) the t-orientation convention
  (toward vs away from the reference point) flips EVERY directional
  component's sign (ratios −1.0000 vs +1.0000 against the reference);
  (d) the moment-arm convention (mid−node vs node−mid) flips every
  moment's sign independently of (c). All four measured
  independently; the reference implementation exhibits (b) itself
  without its own validation keyword.

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
| (A) | Honest stopping criteria in the iterative solver (residual-based; no false-success exits) | **DONE** (lote A+B, 2026-08-28) — success only on the real residual test; breakdown/stagnation false-success exits removed; honest failure (RC≠0) with the real relative residual; `control_solver_bicg_stop -no` as the documented continue-escape. Key result: the 3D flat-residual cantilever (A·b≈0 family) now converges (final error 4.2e-37) — the flat-residual identity is NOT stagnation, CG just needs to keep iterating. |
| (B) | Conjugate Gradient for the SPD system (replacing/augmenting Bi-CG) | **DONE** (lote A+B, 2026-08-28) — runtime symmetry check dispatches CG (symmetric) / honest Bi-CG (non-symmetric: beam dtime-asymmetry, plastic-slip interface with dᵀAd<0); primal-residual monitor (the old monitor measured the transpose residual, which never vanishes on non-symmetric systems — a measurement bug of the dishonest criteria family). |
| (C) | Real monolithic mixed u-sigma solve (MINRES / SuperLU with pivoting) | **NOT NEEDED** (2026-08-28) — the staggered scheme with the element-consistent feedback (D) converges to the element solution in one pass; the monolithic refactor would add no physics for the linear-elastic case. |
| (D) | Element-consistent staggered scheme: the momentum feedback carries the ELEMENT internal force (Bᵀσ_old + dt·K_elem·v) and the σ dofs are recovered by the Lagrange extrapolation of the Gauss-point values (the "same B at the node") | **DONE** (lote C/D, 2026-08-28) — fixed point BEFORE: v* = (dt·K_full)⁻¹·(P − Bᵀσ_old) (the full-constitutive/locked solution, the SRI cancelled); AFTER: v* = (dt·K_elem)⁻¹·(P − Bᵀσ_old) (the ELEMENT solution). Measured: SRI quad4 clamp moment 0.0750 = 0.9375·P·L stable at 32 iterations; plain quad4/quad9/hex8 byte-identical; gforce10/13 converge with the axial N EXACT (1.0000× vs Professional); the arness divergences were input BC bugs (rigid-rotation mechanisms), fixed. Not fixed: the plain quad4/hex8 lock (element physics, opt-in SRI — the hex8 SRI now recovers 0.897× of the EB deflection, 2026-08-29) and the section shear pollution of the raw Q4 σ_xy (pre-existing; the Professional's nodal σ_xy is equally polluted — its exact statics do not come from the raw nodal stress). |
| (E) | Equilibrated sigma state at the end of each step: the ELEMENT_DOF 3D write fixed (the constitutive stress lands at the nder-correct slot) + the LOT-5 2D section kinematics corrected (the inverse of the true Jacobian) | **DONE** (lote 6, 2026-08-29) — gforce7/gforce7_ref N/V/M from 5.0000×/1.2500× to 1.0000×/0.9984×; the 3D reads the real ELEMENT_DOF (the L4 fallback retired); the assembled internal forces equal the loads to the solver tolerance (Bᵀσ = P); suite 209/209. Measured: the "non-equilibrated σ" of the arness was a post-processing artifact — the scheme's fixed point was equilibrated all along. |
| — | Q4/Q8 (hex8) selective reduced integration (opt-in) | DONE — `group_element_selective_reduced_integration` (quad4 2026-08-28; hex8 2026-08-29). SRI hex8 measured: loaded 8×1×1 cantilever 0.897× of the EB deflection (vs 0.221× locked) — close to but below the 2D 0.9375×. Measured caveat for the LINE 2 sensitivity table: the shear-only SRI hex8 has 9 zero eigenvalues isolated (6 rigid + 3 twist) and section-warping zero modes in a mesh (the clamped cantilever system is singular; the load-orthogonal solve still converges) — the classic instability that moved the brick literature to B-bar/ANS; the 2D quad4 SRI is stable. |

## 5. Open questions

- ~~Is the 0.2315× fixed point a property of the staggered scheme itself, or
  an artifact of the dishonest inner solver polluting the outer
  iterations?~~ **ANSWERED (lote A+B)**: it is the scheme's. With the honest
  solver the sweep 1..32 gives identical fixed points (plain 0.0185184 =
  0.23148×, SRI 0.0185298 = 0.23162×) — the 2D inner solve was already
  converging honestly. The road model confirms: with converged solves,
  |M_end| + M_center = 0.143·pL²/8, still impossible statics.
- The hex27 divergence (2.16e+13 after 720 iterations) was not reproduced with
  regular cantilevers; the exact failing mesh from the earlier session is not
  in the repository. (The hex27 ×8 regular cantilever now converges with the
  honest solver: 79 iterations, mom ≈ 1.1·P·(L−x).)
- Moment extraction for meshes with 2+ elements in thickness was
  inconclusive (open question).
- ~~Does the Professional's exact section statics come from the raw nodal
  stress or from something else (DIAG §12.4, question (b))?~~
  **CLOSED (lote 4, 2026-08-28)**: it does NOT come from the raw
  stress AT ALL — not the nodal AND NOT the element integration-point
  stress. Measured in the arness: for the quad9 the recovered nodal
  stresses are the EXACT element-IP averages (the Lobatto recovery is
  the Kronecker delta at the shared nodes), so switching the section
  source from NODE_DOF to ELEMENT_DOF (calcul_force.cc LOT 4) leaves
  the N/V pollution unchanged (gforce7 N 1.24×, V 2.7× — the
  pollution lives in the σ FIELD itself, which the Professional's
  nodal σ shares, DIAG §12.4: clamp σ_xy = −90.9 vs τ = 0). The
  Professional's exact statics (1e-10) are consistent with an
  EQUILIBRIUM/INTERNAL-FORCE-based section calculation ("the element
  forces needed for this option are setup in a timestep", manual
  6.913) — the GNU's section integration over any raw σ field cannot
  reach that exactness. Paper line: "stress-field integration vs
  equilibrium section forces in mixed u-σ formulations".
  **RESOLVED (lote 5, 2026-08-28)**: the equilibrium section forces
  are now IMPLEMENTED in the GNU (calcul_force.cc LOT 5): the section
  resultants = the sums of the element internal forces f_elem =
  ∫Bᵀσ dV of the face nodes = the free-body statics of the loads
  (EXACT for the equilibrated σ states). Measured evidence that this
  IS the Professional's family: (a) the tunnel ring msf_tunnel3d now
  reproduces the Professional's node_dof_calcul DIGIT FOR DIGIT (nor
  0.0998068, shes 0.013727, mom1 −2.0e-4 — the FE-discretized values,
  not the analytic 0.1 of the prescribed field that the stress-field
  integration returned); (b) the fixed-fixed beam identity
  |M_end|+|M_center| = pL²/8 holds at 0.9996 (the GNU's σ-field
  integration gave 0.197×); (c) the 2-element quad4 cantilever gives
  N/V EXACT and M 99.84%. The paper line now has the full contrast:
  stress-field integration (σ-field pollution, formulation-dependent)
  vs equilibrium section forces (free-body statics, formulation-
  independent — even the shear-locked quad4's equilibrated solution
  yields the exact statics). LIMITATION measured: the equilibrium
  resultants are the free-body statics ONLY for equilibrated σ states;
  the GNU's coarse multi-step quad9/hex8 runs carry a NON-equilibrated
  recovered σ (pre-existing staggered-scheme property — the L4
  baseline measured it as N 1.24×/V 2.7× by field integration; the
  equilibrium resultants amplify it to 5×/1.25×). The Professional
  produces the equilibrated σ (its solver), which is why its statics
  are exact on the same meshes — closing THAT gap is a SOLVER issue
  (making Bᵀσ = P hold), not a section-post-processing one.

---

## 6. The definitive comparison: Tochnog Professional (2026-08-28)

The user obtained the Tochnog Professional binary (version 02-08-2026, from
the author's public Drive "PublicDennis") and we ran it on the SAME
pathological models. **The bug does not survive in the Professional.**

- **Cantilever quad4 1-element-in-thickness** (the Professional's own
  `force7.dat` statics rebuilt with quad4 instead of quad9): Professional
  gives N = -12.34, V = +100, M(x=50) = -5000, M(x=0) = -10000,
  M(x=100) = 0 — ALL EXACT to 1e-10. The GNU gives 0.2315×.
- **Fixed-fixed beam, uniform load, quad4 1-in-thickness** (the user's
  road-base case): |M_end| + |M_center| = 825 + 425 = **1250 = pL²/8
  EXACT** to 1e-10. The individual 825/425 are the correct deep-beam
  Timoshenko fixed-end moments for h/L = 0.1 (Euler-Bernoulli would give
  833/417); the sum is the statics identity the GNU violates by ≈7×.
- The Professional's own validation family uses quad9 (force7/8) and 2×2
  in-section hex8 (force10/13) — defensive test design; the quad4
  1-in-thickness case itself works correctly in the Professional.
- The GNU's staggered "stress follows the principal unknowns" architecture
  is documented in the 2011 open manual — the defect is inherited
  open-source lineage, not present in the proprietary product.

Full evidence: `ProjectDocs/DIAG-SOLVE-MIXTO.md` §11.

## 7. Status of the three lines (updated)

| Line | Status | New material from the Professional comparison |
|------|--------|-----------------------------------------------|
| 1 — staggered fixed point + solver honesty | Evidence complete; repairs A+B, C/D and E done (2026-08-28/29) | The Professional proves the correct answer is attainable with the same element/mesh — isolates the open-source scheme as the defect carrier; repair E closes the arness gap (Bᵀσ = P at the end of every step) |
| 2 — solver formulation sensitivity | Initial measurements; needs the repair to complete the comparison table | Bi-CG ≡ SuperLU (md5) + the Professional's exact results bracket "what the solver should deliver" |
| 3 — pedagogical hand calculation | Case study complete | The road-base case now has a positive control: the SAME model gives exact statics in the Professional — the detector works, the defect is software-specific |

