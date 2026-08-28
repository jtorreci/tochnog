# group_element_selective_reduced_integration

## The phenomenon: shear locking of the bilinear quad4 (parasitic shear)

This is a CLASSIC result of the finite element literature (T.J.R. Hughes,
"The Finite Element Method: Linear Static and Dynamic Finite Element
Analysis", chapter on reduced/selective integration; see also "parasitic
shear" of the 4-node bilinear element in bending). It is documented here
because it was first *measured* in this repository during the
`materi_stress_force` sprint (lote 2, 2026-08-28), and because the opt-in
fix implemented with this keyword is only understandable with the
mechanism in mind.

### Mechanism

A cantilever beam under bending has a quadratic transverse displacement
field. For a beam along x with curvature κ:

```
u_y(x) ~ -κ·x²/2        (quadratic)
u_x(x,y) = -κ·x·y        (bilinear — exactly representable)
```

The bilinear shape functions of the quad4 CANNOT represent the quadratic
`u_y`. When the element is forced into bending, the best bilinear fit
cannot make the shear strain vanish:

```
γ_xy = ∂u_x/∂y + ∂u_y/∂x = -κ·x + (linear fit slope of u_y)  ≠ 0
```

The element responds by generating a SPURIOUS shear strain γ_xy ≠ 0
(parasitic shear). The constitutive law converts it into spurious shear
energy `∫G·γ² dV`, so the bending response of the element becomes TOO
STIFF: the element under-estimates moments and displacements. The effect
is independent of the mesh refinement along the beam (it is a property of
the single element in the thickness direction), which is why it is called
*locking*.

### Measured in this repository (lote MSF L2, 2026-08-28)

Cantilever 2D, plane stress, 8 quad4 elements along the length, ONE
element in the thickness, L=8, h=1, E=1000, ν=0.3, tip load P=1e-2
(test `qsri_beam2d`):

| quantity | quad4 (full integration) | analytic | ratio |
|----------|--------------------------|----------|-------|
| section moment at the clamp | 0.0185 | P·8 = 0.08 | **0.231×** |
| section shear (interior) | 0.00741 | P = 0.01 | **0.741×** |

The same model with quad9 gives the section moment within 1-3% of P·(8−x)
(`msf_beam2d`). The same phenomenon exists for the trilinear hex8 in 3D
(`msf_sheet3d_hex8`, lote MSF L3): the hex8 cannot represent the
quadratic bending field `u_x ~ y²` either; the hex8 SRI is FUTURE WORK
(only the 2D quad4 is wired to this keyword).

## The fix: selective reduced integration (SRI, Hughes)

The classic fix removes the parasitic shear from the stiffness without
introducing spurious zero-energy modes:

- the NORMAL/volumetric terms of the constitutive matrix
  (`σ_xx`, `σ_yy` and the ν-coupling, i.e. everything that is NOT shear)
  are integrated with the FULL rule,
- the SHEAR term (γ_xy) is integrated with ONE Gauss point at the element
  centroid.

For the pure-bending mode, the parasitic shear strain vanishes AT THE
CENTROID of the element (the optimal bilinear fit of the quadratic `u_y`
has exactly the right slope at the centroid), so the spurious shear
energy disappears from the stiffness, while the normal terms (integrated
fully) keep the element stable: no hourglass mode appears.

### IMPORTANT codebase finding (measured during this work)

The classic SRI result (moment ≈ P·(8−x)) assumes the "full" rule to be
the standard 2×2 GAUSS rule. The default quadrature of this GNU for the
quad4 is the 2×2 LOBATTO rule at the element CORNERS (polynom.cc:
`integration_method=-LOBATTO` with `integration_points=-MAXIMAL`). The
corner rule OVER-integrates the quadratic bending energy (1.5×), which
would cap the benefit of a shear-only SRI at ≈ 0.35× of the exact moment
(measured). Therefore the keyword also switches the full rule of the
group to the classic 2×2 GAUSS rule. Both decisions are documented in
the measured values below.

## Keyword

```
group_element_selective_reduced_integration <element_group> -yes
```

- Opt-in, DEFAULT OFF: without the keyword (or with `-no`) the group
  behaves EXACTLY as before (byte-identical path; the full 191-run suite
  stays green unchanged).
- Scope (validated with the keyword): 2D bilinear **quad4** (nnol==4),
  LINEAR ELASTICITY only, no axisymmetry, no `materi_displacement`
  (large displacement), no plasticity/damage/maxwell models. Any other
  use is ignored with a one-time warning (see `sri_quad4_active()`).
  hex8 SRI in 3D is future work.

## Files and functions

- `miscel.cc` — `sri_quad4_active(element, element_group, name, nnol)`:
  the SINGLE scope gate used by BOTH consumers below (so the quadrature
  and the stiffness split can never disagree). Reads the keyword
  (`GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION`) and applies the
  exclusions above.
- `polynom.cc` — `pol()`: when `sri_quad4_active()` is true, the
  integration method of the quad4 is forced to `-GAUSS` (the classic
  2×2 full rule for the normal terms; the codebase default is the 2×2
  Lobatto corner rule).
- `materi.cc` — `materi()`: the constitutive split. At every
  integration point, after `set_stress()` fills `ddsdde_total`:
  ```c
  if ( sri_quad4 ) {
    sri_g = ddsdde_total[1*MSTRAIN+1];   // G: gamma_xy entry (stress_indx(0,1))
    ddsdde_total[1*MSTRAIN+1] = 0.;      // D_norm for the full integration
  }
  matrix_atba( new_b, ddsdde_total, stiffness, work, MSTRAIN, nnol*ndim );
  ```
  and the reduced shear stiffness is assembled at the element centroid:
  the bilinear local derivatives at iso (0,0) are derived from the
  index convention of polynom.cc (row-major from the bottom: node
  `i` at local `((2*(i%2)-1), (2*(i/2)-1))`, NOT the textbook
  counter-clockwise order), the Jacobian at the centroid handles
  distorted elements, and the 1×1 Gauss weight is 4:
  ```c
  stiffness_shear[i*nnol*ndim+j] = sri_g * sri_volfac *
    4.*detj_center * b_shear[i] * b_shear[j];   // b_shear = gamma_xy row of B
  ```
  The reduced term is added to `element_matrix` scaled by `1/npoint`
  (materi() is called once per integration point; the npoint calls sum
  exactly to the full reduced integral). The `element_lhside` diagonal
  gets the same treatment.
- `database.cc` — keyword registration (type `INTEGER`, length 1,
  `data_class = GROUP_TYPE`, after `GROUP_INTEGRATION_POINTS`).
- `check.cc` — requires `materi_stress` and `materi_velocity`.
- Enums `GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION` in `tochnog.h` /
  `tochnog-mod.h` (in sync).

## Verification (family `qsri`, all A/B, see build_safe.sh)

- `qsri_beam2d` / `qsri_beam2d_sri` (materi_stress_force.440/441):
  the cantilever 8×1 above. WITHOUT SRI the documented lock is
  reproduced EXACTLY (mom = 0.231×, she = 0.741×). WITH SRI: mom =
  0.312× and she = 0.65× at the same sections — a real improvement,
  robust across the Bi-CG and the direct solver.
- `qsri_patch_s` / `_off` (442/443): simple shear, PRESCRIBED linear
  field γ = 1e-3, single quad4: she = G·γ = 0.384615384615 EXACT
  (12 digits), nor = mom = 0, byte-identical values with and without
  SRI — the SRI integrates constant states exactly.
- `qsri_patch_t` / `_off` (444/445): uniaxial tension, prescribed
  ε_xx = 1e-3: nor identical with and without SRI (1.025641; the small
  deviation from the hand value E'·ε·h is the σ-assembly of the mixed
  scheme, identical in both runs).
- `qsri_modes_rigid` (449): prescribed rigid ROTATION → stresses
  EXACTLY zero (the 3 rigid modes — 2 translations + 1 rotation — are
  exact zero modes of the SRI stiffness).
- `qsri_modes` (448): prescribed hourglass deformation u_y = a·x·y →
  NONZERO stresses (the 4th mode has positive energy: no zero-energy
  hourglass mode).
- Eigenvalue analysis of the isolated SRI quad4 (documented, computed
  independently): exactly 3 zero eigenvalues (rigid modes), the 4th
  mode positive, rank 5/8. The classic 1-point (total) reduced
  integration would leave the hourglass mode at zero energy; the
  selective scheme does not.

## Solver note (GOTCHA of the GNU mixed scheme, documented)

The classic SRI result (moment ≈ P·(8−x) for the 8×1 cantilever) is
recovered when the u-system `K_SRI·u = P` is solved reliably: the exact
solution of the C-assembled system (verified with Gaussian elimination on
the extracted element matrices) gives tip deflection 99.3% and clamp
moment 93.75% of the analytic values. The GNU solve of this mixed u-σ
formulation, however, is numerically unreliable for the 1-element-in-
thickness quad4 system (condition number ≈ 2·10⁴): the Bi-CG stops at a
misleading residual and the direct solvers (which renumber the nodes)
produce null-mode-polluted fields (measured). This is the SAME family of
GOTCHA as the 3D mixed-solve degeneration documented in the
materi_stress_force sprint — it is a solver property, NOT a defect of
the SRI element (whose stiffness is verified exact). The measured A/B
above (0.231× → 0.312×) is therefore the reproducible improvement in the
current solver, and the classic 1.0× value is documented as the
solver-independent reference.

## Hardcoded parameters / pending refactorings

- The shear term index is hardcoded as `1*MSTRAIN+1` (stress_indx(0,1),
  the γ_xy entry); correct for the 2D plane formulation (isotropic: the
  shear does not couple with the normal rows).
- The reduced point is the centroid of the isoparametric element (1×1
  Gauss, weight 4) — valid for any distorted quad4 (the Jacobian at the
  centroid is used).
- The centroid B is built from the reference coordinates passed to
  materi() (== the deformed ones when `materi_displacement` is off,
  which is a scope condition).
- hex8 SRI in 3D: future work (same mechanism, the shear terms are
  γ_xy, γ_xz, γ_yz).
