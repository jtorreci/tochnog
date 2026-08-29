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
quadratic bending field `u_x ~ y²` either (measured in this extension:
0.221× of the Euler-Bernoulli deflection for the loaded 8×1×1
cantilever without the keyword). The hex8 SRI is implemented in this
extension (2026-08-29, section below).

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
- Scope (validated with the keyword): 2D bilinear **quad4** (nnol==4)
  and 3D trilinear **hex8** (nnol==8), LINEAR ELASTICITY only, no
  axisymmetry, no `materi_displacement` (large displacement), no
  plasticity/damage/maxwell models. Any other use is ignored with a
  one-time warning (see `sri_active()`).

## Files and functions

- `miscel.cc` — `sri_active(element, element_group, name, nnol)`
  (renamed from `sri_quad4_active` with the hex8 extension): the SINGLE
  scope gate used by BOTH consumers below (so the quadrature and the
  stiffness split can never disagree). Reads the keyword
  (`GROUP_ELEMENT_SELECTIVE_REDUCED_INTEGRATION`) and applies the
  exclusions above.
- `polynom.cc` — `pol()`: when `sri_active()` is true, the
  integration method of the element is forced to `-GAUSS` (the classic
  2×2 / 2×2×2 full rule for the normal terms; the codebase default is
  the Lobatto corner rule).
- `materi.cc` — `materi()`: the constitutive split. At every
  integration point, after `set_stress()` fills `ddsdde_total`:
  ```c
  if ( sri_on ) {
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
misleading residual. This is the SAME family of
GOTCHA as the 3D mixed-solve degeneration documented in the
materi_stress_force sprint. The measured A/B
above (0.231× → 0.312×) is therefore the reproducible improvement in the
current solver, and the classic 1.0× value is documented as the
solver-independent reference.

**Correction (2026-08-28, DIAG-SOLVE-MIXTO.md)**: the diagnosis was
refined — the linear solver is NOT the limiter. The Bi-CG system is
velocity-only and SPD (the σ dofs never enter the global matrix; the σ
field is advanced by a lumped staggered update), and a direct solve
(`control_options_solver -matrix_superlu`) of the SAME assembled system
is **byte-identical** to Bi-CG (md5 `62e82a92...`, `mom = 0.0249998465889`
at the default 2 equilibrium iterations). The true limiters are: (a) the
Bi-CG exit criteria, which report *success* on breakdown/stagnation with
`error ≫ check_error` (wrong solutions with `RC=0`), and (b) the
staggered fixed point, which converges to `mom = 0.2315×` for BOTH the
plain and the SRI quad4 (the `0.3125×` SRI value is a transient of the
2-iteration default; the σ-gradient RHS re-introduces the full
constitutive shear at equilibrium, cancelling the SRI). The `93.75%`
reference is the pure-displacement system `K_SRI·u = P`, which the GNU
mixed scheme does not assemble.

## HEX8 extension (3D, implemented 2026-08-29) — same keyword, same scope

The keyword now applies to the trilinear **hex8** (3D) as well, with the
same opt-in default-OFF scope: linear elasticity, no
`materi_displacement`, no plasticity/damage/maxwell (axisymmetry is a 2D
concept and does not apply). The single gate is `sri_active()`
(miscel.cc, renamed from `sri_quad4_active()`): `(QUAD4 && nnol==4 &&
ndim==2) || (HEX8 && nnol==8 && ndim==3)`.

### The 3D split

`D = D_norm + D_shear` with the shear diagonal entries of the 6×6
isotropic tangent zeroed for the full rule: indices `stress_indx(0,1)=1`
(γ_xy), `stress_indx(0,2)=2` (γ_xz) and `stress_indx(1,2)=4` (γ_yz).
The reduced shear stiffness is the 1×1×1 Gauss point at the centroid
(weight 8), built from the three engineering shear rows of B at the
centroid (the polynom.cc:549-599 convention, local derivatives derived
from the index, the transpose-Jacobian physical derivatives):

```
gamma_xy: [dN/dy, dN/dx, 0]      gamma_xz: [dN/dz, 0, dN/dx]
gamma_yz: [0, dN/dz, dN/dy]
K_shear = volfac * 8*detJ_c * ( G_xy*b_xy*b_xy^T + G_xz*b_xz*b_xz^T
                                + G_yz*b_yz*b_yz^T )
```

### The element-consistent feedback (D-c) in 3D

The momentum RHS must use the SAME reduced shear as the matrix (fix D-c
of DIAG-SOLVE-MIXTO §12, without which the SRI is cancelled at the
fixed point — measured for the quad4). For the hex8 the current-iterate
shear increments of the feedback stress are zeroed in `sigvec` and the
reduced internal force `-dt*K_shear*v` is added to the RHS, exactly as
for the quad4. **Indexing gotcha (measured)**: `inc_ept` is the
MDIM×MDIM TENSOR strain (`inc_ept[idim*MDIM+jdim]`) while `sigvec` is
Voigt-indexed (`stress_indx`). In 2D the indices coincide; in 3D the yz
entry differs (tensor `1*MDIM+2 = 5` vs Voigt `stress_indx(1,2) = 4`),
so the strain MUST use the tensor index `idim*MDIM+jdim`.

### The consistent σ recovery (D-b) in 3D

The normal stress components (`stress_indx` 0, 3, 5 — the same list as
the 2D) get the trilinear Lagrange extrapolation of the 2×2×2
Gauss-point values to the nodes (`sri_stress_recovery_weight` with the
3D tensor product `w_xi*w_eta*w_zeta`, node iso coordinates
`(2*(i%2)-1, 2*((i/2)%2)-1, 2*(i/4)-1)`); the shear components keep the
h-weighting. The section (calcul_force.cc) switches the full rule to
Gauss 2×2×2 (`msf_element_rule`) and replaces the full-rule shear part
of the element internal forces by the reduced 1-point shear internal
force (the mean IP shear stresses, the moduli cancel — the same
correction as the 2D quad4, extended to the three shear rows).

### VERIFICATION — a measured deviation from the classic claim

**The shear-only SRI hex8 is NOT hourglass-free** (unlike the 2D
quad4). Verified with the eigenvalue analysis of the isolated element
(independent reconstruction of B and D in the code conventions):

| rule | zero eigenvalues of the isolated hex8 |
|------|---------------------------------------|
| full 2×2×2 (everything) | 6 (the rigid modes) |
| **SRI (2×2×2 normals + 1-point shear)** | **9 (6 rigid + 3 twist)** |
| 1-point (fully reduced) | 18 (6 rigid + 12 hourglass) |

The three twist modes are `u = (y·z,0,0)`, `u = (0,x·z,0)`,
`u = (0,0,x·y)`: all normal strains zero, the shear strains vanish at
the section centroid, so the 1-point rule misses them. In a mesh the
additional **section-warping** modes appear, `u_y = A(x)·(2z−1)`,
`u_z = A(x)·(2y−1)` with A a piecewise-linear hat (measured: 8 such
modes in the clamped 8×1×1 cantilever — the clamped system is
SINGULAR). The 2D quad4 SRI is stable because its warping analog
`u_y = A(x)·(2y−1)` has a normal strain `ε_yy = 2A(x)` that the full
normal rule stabilizes; the 3D warping modes have all normal strains
zero. This is the KNOWN reason the FE literature moved from the
shear-only selective integration of the 8-node brick to the B-bar /
assumed-strain (ANS) formulations.

**The loaded cantilever still works in the GNU mixed scheme** (measured,
family `qsri3d`): the 8×1×1 hex8 cantilever (L=8, section 1×1,
P=1e-2 tip, E=1000, ν=0.3 — the direct 3D analog of `qsri_beam2d`):

| quantity | hex8 full (OFF) | hex8 SRI (ON) | Euler-Bernoulli |
|----------|-----------------|---------------|-----------------|
| tip deflection | 0.00453 (0.221×) | 0.01838 (**0.897×**) | 0.02048 |
| axial σ_zz at the clamp | ±0.123 (0.26×) | ±0.509 (1.06×) | ±0.48 |
| section mom1s (fallback σ-field) | 0.0788 (0.985×) | 0.0909 (1.14×) | P·L = 0.08 |

The OFF value reproduces the documented hex8 lock (the 0.23× family of
MSF L3); the SRI recovers ~90% of the Euler-Bernoulli deflection and
~100% of the bending stress — close to but below the 2D quad4 SRI's
93.75% (the task-expected "the hex8 does not reach the exact 2D
0.9375"; measured: 0.897×). The solve converges because the load is
orthogonal to the warping nullspace and the honest CG stays out of it;
meshes whose BCs do NOT kill the warping modes (free lateral faces,
simply-supported configurations) are at risk — the documented
limitation of the shear-only SRI hex8.

Patch tests (constant states, single/few elements, prescribed
velocities): 3D simple shear `she = G·γ·t = 0.5` EXACT and 3D confined
tension `nor = 1.346154` EXACT, with IDENTICAL values with and without
the SRI (the SRI integrates constants exactly) — `qsri3d_patch_s/_off`
(700-703) and `qsri3d_patch_t/_off` (704-707). Rigid rotation: section
0 within 1e-17 (the 6 rigid modes are exact zeros of the SRI matrix) —
`qsri3d_modes_rigid` (708). Prescribed twist `u_x = a·(y−0.5)·(z−0.5)`:
the SRI section reads ~0 (the zero-energy mode is invisible to the
reduced rule) while the full constitutive law still sees the shear at
the 2×2×2 points (the nodal σ output is nonzero) — the DOCUMENTATION
test of the twist mode, `qsri3d_modes` (709). Cantilever A/B:
`qsri3d_beam` (710, SRI) / `qsri3d_beam_off` (711) with the tip
deflection targets in the .dat (0.18379 vs 0.04528, discriminating).

**Arness (Professional comparison)**: the hex8 cases gforce10/gforce13
(default, no SRI — they use `materi_displacement`, outside the SRI
scope anyway) stay at the axial N = 12.34 EXACT (1.0000× vs the
Professional) — no regression. **Suite**: 209/209 runs + all file
checks green on a clean build (201 pre-existing + 8 new); the default
(no keyword) path is byte-identical (the SRI gates are all
`GET_IF_EXISTS` + element-scope checks).

## Hardcoded parameters / pending refactorings

- The quad4 shear term index is hardcoded as `1*MSTRAIN+1`
  (stress_indx(0,1), the γ_xy entry); the hex8 shear diagonal indices
  are hardcoded as `1*MSTRAIN+1`, `2*MSTRAIN+2`, `4*MSTRAIN+4`
  (isotropic: the shear rows of D are decoupled, so the diagonal
  zeroing is the exact D_norm; a general anisotropic tangent would
  need the full row/column zeroing).
- The 3D strain indexing gotcha: `inc_ept` is the MDIM×MDIM tensor
  strain, so the yz shear entry is `inc_ept[1*MDIM+2]` (NOT
  `inc_ept[stress_indx(1,2)] = inc_ept[4]`) — the quad4 2D coincidence
  does not hold in 3D.
- The reduced point is the centroid of the isoparametric element
  (1×1 / 1×1×1 Gauss, weight 4 / 8) — valid for any distorted element
  (the centroid Jacobian is used).
- The centroid B is built from the reference coordinates passed to
  materi() (== the deformed ones when `materi_displacement` is off,
  which is a scope condition).
- KNOWN LIMITATION (measured 2026-08-29): the shear-only SRI hex8
  retains zero-energy modes (3 twist modes isolated, section-warping
  modes in a mesh — the clamped cantilever system is singular). The
  load-orthogonal cantilever still solves (0.897× EB), but general
  meshes need the stabilization that the literature moved to
  (B-bar / assumed-strain). The 2D quad4 SRI does NOT have this
  limitation (its warping mode has a normal strain).
- The section (calcul_force.cc) reads the element internal forces of
  the ELEMENT_DOF σ field with the documented fallback to the
  recovered nodal σ when the staggered loop does not propagate the
  strains (hex8 + `derivatives`); the fallback σ of the loaded hex8 is
  NOT in equilibrium, so the beam section values (710/711) are the
  polluted σ-field readings — the deflection targets are the
  discriminating checks.



- The shear term index is hardcoded as `1*MSTRAIN+1` (stress_indx(0,1),
  the γ_xy entry); correct for the 2D plane formulation (isotropic: the
  shear does not couple with the normal rows).
- The reduced point is the centroid of the isoparametric element (1×1
  Gauss, weight 4) — valid for any distorted quad4 (the Jacobian at the
  centroid is used).
- The centroid B is built from the reference coordinates passed to
  materi() (== the deformed ones when `materi_displacement` is off,
  which is a scope condition).
- hex8 SRI in 3D: IMPLEMENTED (2026-08-29, the section above; the
  shear terms γ_xy, γ_xz, γ_yz at the 1×1×1 centroid point) — with the
  measured zero-energy twist/warping modes documented as the known
  limitation of the shear-only selective integration of the brick.
