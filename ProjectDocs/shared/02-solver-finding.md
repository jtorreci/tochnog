# The solver finding: a 30-year-old defect in the open-source line

## In one paragraph

The Tochnog GNU line solves the momentum equation with stress as a
non-principal unknown that "follows" the velocities through a **staggered
(operator-split) scheme**. We found that this scheme's **fixed point is not
the displacement solution**: for a 1-element-in-thickness bending beam it
converges to section moments of **0.23× the exact statics** — physically
impossible results that pass silently because the Bi-CG solver's **stopping
criteria declare success without converging**. The defect is present in the
open-source line since at least the 2009-era architecture; it does **not**
exist in Tochnog Professional (the same element, the same mesh, gives exact
statics there). We diagnosed it with instrumented evidence, fixed it, and
validated the fix against the Professional binary.

## The symptom (how it was found)

While implementing `post_calcul -materi_stress -force` (section forces and
moments), a real-world check surfaced: in a 2D plane-strain road base model
(subgrade + concrete slab), integrating the concrete stresses gave
|M_end| + |M_center| ≈ 0.23·pL²/8, far below pL²/8 — a violation of exact
statics (for a fixed-fixed beam |M_end| + |M_center| = pL²/8 is exact,
independent of the material). The quad9 element gives the correct answer;
the quad4 does not — the classic bilinear shear locking, but amplified by
the scheme.

## The mechanism (measured, with code references)

1. **The "mixed" system is not mixed**: the stress dofs never enter the
   global matrix (`input.cc`: no `dof_principal` for stress). The solved
   matrix is velocity-only, symmetric positive definite. Stress advances
   through a diagonal lumped update and its `−Bᵀσ` term returns to the
   velocity right-hand side.
2. **The fixed point is wrong**: iterating the outer scheme to convergence
   gives `v* = (dt·K_full)⁻¹(P − Bᵀσ_old)` — the **full-constitutive
   (locked) solution, always**. The element matrix (even with selective
   reduced integration) only shapes the transient; it cancels at the fixed
   point. Measured: section moment converges to 0.2315× for both plain
   quad4 and SRI quad4.
3. **The Bi-CG stopping criteria are dishonest**: three of the four exits
   report success without convergence — a breakdown `|dAd| < 1e-16`
   returns `x = 0`; a stagnation test returns `RC = 0` with error 4.7e-7
   against a tolerance of 1e-12. For symmetric structures the first CG step
   lands exactly on a flat residual plane (identity `|r|²|Ād|² = 2(dᵀĀd)²`
   verified to 17 digits), which the old criteria mistook for convergence.
4. **The linear solver was never the problem**: SuperLU and Bi-CG produce
   byte-identical results (same md5). Changing solver changes nothing — the
   scheme is the defect carrier.

## The fixes (all committed, 201-test suite green)

- **(A) Honest stopping criteria**: success is reported only on the real
  residual test; breakdown and stagnation no longer fake success; failures
  are honest (`RC ≠ 0` with the true relative residual).
- **(B) CG for the SPD system**: a runtime symmetry check dispatches to
  plain CG for symmetric systems (velocity/temperature/pressure) and to the
  honest Bi-CG for genuinely non-symmetric ones (beam coupling,
  plastic-slip interfaces).
- **(C/D) Consistent fixed point**: the shear feedback of the staggered
  scheme is integrated like the element matrix, and the stress recovery
  uses a consistent (superconvergent) extrapolation for the normal
  components. The fixed point is now **the element's own solution** —
  convergence in one pass, no drift.
- **(SRI) Quad4 selective reduced integration** (opt-in keyword
  `group_element_selective_reduced_integration`): the classical Hughes fix
  for the bilinear locking, verified at element level (exact stiffness,
  3 rigid modes, no hourglass, exact patch tests). End-to-end it now gives
  the classic reference: **mom = 0.9375·P·L** (the 1-element-in-thickness
  cantilever), stable at 32 iterations.
- **(L5) Section statics by internal forces**: the section-force
  post-processor computes forces from the element internal force vector
  (free-body equilibrium), not from the polluted stress field. Result:
  exact statics regardless of the element formulation — even the locked
  quad4 gives the exact section forces (the lock remains in the deflection
  and the stress field, not in the section resultants).

## Validation against Tochnog Professional

| Model | Quantity | GNU (before) | GNU (after) | Professional |
|---|---|---|---|---|
| Cantilever quad4, 1-in-thickness | M(x=L/2) | 0.23× (locked fixed point) | 0.9984× (equilibrium) | 1.0000× (exact) |
| Fixed-fixed beam, uniform load | \|M_e\|+\|M_c\| / pL²/8 | 0.197× (impossible) | 0.9996× | 1.0000× (exact) |
| Tunnel ring (hex27) | nor, shes, mom1 | — | **digit-for-digit identical** to the Professional | reference |
| Cantilever quad9 | M | 0.995× | 0.998× | exact |

The Professional's own validation tests use `post_calcul -materi_stress
-force` with exact-statics targets; our implementation reproduces them.

## Why the suite never caught it

The test suite is a regression/coverage suite: targets were calibrated from
the code's own output, and no legacy test checks an integrated section
moment in a 1-element-in-thickness bending configuration. Wrong answers
that "look reasonable" (0.23×, not NaN) pass forever. A one-line hand
calculation — the statics check |M_e| + |M_center| = pL²/8 — is what
detects it. (This is the subject of one of our paper lines; see
`03-papers.md`.)

## Full documentation

The paper-grade diagnosis (with all evidence: eigenvalue spectra, matrix
dumps, iteration sweeps, code references) is in
`ProjectDocs/DIAG-SOLVE-MIXTO.md` (English, 600+ lines). The comparison
harness and tables are in `ProjectDocs/VALIDACION-PROFESIONAL.md`.
