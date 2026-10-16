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
group_element_selective_reduced_integration <element_group> -yes|-volumetric|-lambda
```

GNU-only extension (no Professional counterpart). `-yes` = the shear
terms at 1 point (the original feature, unchanged); `-volumetric` (C80.1,
2026-10-16) = the volumetric part at 1 point (B-bar); `-lambda` (C80.2,
2026-10-16) = the λ·m·mᵀ part at 1 point (Hughes 1977 λ/μ split). See
the sections "C80.1" and "C80.2" below. The C80.1 value `-both` was
removed in C80.2 and is rejected by the input reader.

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

## C80.1 — `-volumetric` (selective volumetric integration, 2026-10-16)

Owner decision 2026-10-16: extend this keyword (no new keyword) with
`-volumetric` (and `-both`, shear + volumetric, removed in C80.2);
opt-in, never the default; GNU-only.

### Formulation (chosen and why)

Volumetric/deviatoric split of the tangent, with the volumetric part at
the centroid (Hughes 1980 selective integration / B-bar with the
centroid dilatation, small strain):

```
D = D_dev + kappa * m m^T,   m = (1,1,1,0,0,0) (Voigt stress_indx 0,3,5)
kappa = m^T D m / 9          (= lambda + 2G/3 for isotropic elasticity)
K = sum_GaussIP B^T D_dev B w detJ  +  kappa * volfac * W * detJ_c * b b^T
b[inol*ndim+idim] = dN_inol/dx_idim at the centroid,  W = 4 (quad4) / 8 (hex8)
```

Why this split: it is the volumetric part the owner specified, it is
the split that carries over to plastic incompressibility (the B-bar
standard for J2 plasticity), and the deviatoric part keeps every
hourglass mode stiff. The alternative, the λ/μ split of Hughes 1977
(reduce only `λ m m^T`), was measured with an independent solve: the
same Lamé accuracy (−0.24 % at ν=0.4999) and a better bending response
on coarse meshes (cantilever 10×1, ν=0.3: −4.9 % vs +22 %). It was NOT
chosen, because it is not a volumetric/deviatoric split and does not
extend to plasticity. The owner then added it as `-lambda` (C80.2,
below) next to `-volumetric`.

Plane strain: ε_zz = 0, so the dilatation is ε_xx+ε_yy, and σ_zz gets
κ·θ_c plus its deviatoric part. Plane stress (membrane) is refused
(volumetric bit dropped, see below). Axisymmetry is refused by the
existing scope gate (the hoop strain u_r/r would enter the centroid
dilatation and the volume factor 2πr; not implemented).

### Files and functions

- `miscel.cc` `sri_active()` now returns a MODE MASK: `SRI_SHEAR` (1)
  for `-yes`, `SRI_VOLUMETRIC` (2) for `-volumetric`, `SRI_LAMBDA` (4)
  for `-lambda` (C80.2), 0 otherwise (macros next to the prototype in `tochnog.h` /
  `tochnog-mod.h`). All boolean consumers (`pol()`, `delete.cc`,
  `general.cc`, the msf rule) are unchanged: every mode switches the
  full rule to Gauss and keeps the SRI recovery weights (Lagrange for
  the normal components, h for the shear). For a membrane group the
  volumetric bit is dropped. A value other than -no/-yes/-volumetric/
  -lambda is a `db_error`.
- `materi.cc`:
  - `sri_centroid_derivatives()` (static): physical dN/dx at iso 0, the
    centroid detJ and coordinates. It uses the same index convention
    and transpose Jacobian as the shear blocks, for both quad4 and hex8.
  - `materi()`, in the split block after the shear split:
    - κ from `ddsdde_total`;
    - θ_ip = tr(inc_ept);
    - θ_c = dt·b·v, from `new_dof` velocities, the same source as the
      shear fix D-c.
    - `sigvec[0,3,5] -= κ·θ_ip`: the D-c feedback. The IP volumetric
      increment leaves the momentum RHS.
    - `new_sig[aa] += κ·(θ_c − θ_ip)`: the STORED stress is the B-bar
      stress, used by the recovery, element_dof and the next step's old
      stress.
    - `ddsdde_total[0,3,5 × 0,3,5] -= κ` (D_dev).
  - The reduced stiffness block adds K_vol to `stiffness_shear`, which
    now carries the whole reduced stiffness and is zeroed first for
    `-volumetric`. The existing matrix (`dt·K/npoint`) and RHS
    (`−dt·K·v/npoint`) loops then treat shear and volume alike.
  - The warning distinguishes `-yes` from `-volumetric`/`-lambda`; it is
    printed when the applied mask differs from the requested one.
- `calcul_force.cc`: the section internal force replaces the full-rule
  shear only for the shear bit (`sri_active() & SRI_SHEAR`). For the
  volumetric part, the stored B-bar stress integrated with the Gauss
  rule IS the reduced internal force, because Σ_ip dN/dx·detJ·w = 4·
  detJ_c·dN/dx(0) exactly for the quad4 (the integrand is bilinear in
  ξ, η) and for the parallelepiped hex8.
- `database.cc`: the names `volumetric` and (C80.2) `lambda`. The enums
  `VOLUMETRIC` and `LAMBDA` are appended before `LAST_DUMMY` in both
  headers (in sync); `LAMBDA` took the slot of the removed `BOTH`, so
  no enum value moved.

### Verification (measured 2026-10-16)

Independent check: a pure-displacement Python solve of the same
element (2×2 Gauss D_dev + 1-point κ b b^T) reproduces the GNU
staggered results to all printed digits on every cantilever (Gauss,
`-yes`, `-volumetric`). For linear elasticity, the D-c
feedback makes the fixed point converge in one pass, as for the shear.
The Gauss, `-yes` and default runs of the study are byte-identical to
the previous binary (51/51 `.dbs`).

Accuracy against analytic solutions (plane strain; Lamé a=1, b=2, p=1;
Timoshenko-Goodier cantilever with E'=E/(1−ν²), ν'=ν/(1−ν), exact
clamp displacements). The table gives the u_r(a) error | the max
deviation of the IP mean stress from p a²/(b²−a²), in % (the
checkerboard measure):

| Lamé 8×16 | Gauss | `-yes` | `-volumetric` | Pro f=2 | Pro f=10 |
|-----------|-------|--------|---------------|---------|----------|
| ν=0.3 | −0.33 / 31 | −0.26 / 32 | −0.17 / 4 | −0.25 / 27 | −0.01 / 5 |
| ν=0.49 | −4.01 / 606 | −3.94 / 607 | −0.19 / 4 | −3.04 / 527 | −0.13 / 105 |
| ν=0.4999 | −80.16 / 12531 | −80.16 / 12533 | −0.19 / 4 | −75.18 / 13493 | −10.80 / 9401 |

The hex8 Lamé (one layer, w=0) is identical to the quad4 to all
digits. The GNU Gauss column equals the Pro `-gauss` column of the C80
study exactly.

Cantilever, tip deflection error:

| case | Gauss | `-yes` | `-volumetric` | λ-split (C80.2 `-lambda`) |
|------|-------|--------|---------------|---------------------------|
| ν=0.3 10×1 | −36.28 | −18.32 | +22.08 | −4.86 |
| ν=0.3 20×2 | −12.42 | −5.28 | +4.70 | −1.25 |
| ν=0.3 40×4 | −3.44 | −1.38 | +1.14 | −0.32 |
| ν=0.49 10×1 | −91.23 | −91.09 | +67.15 | +30.30 |
| ν=0.49 20×2 | −70.61 | −70.24 | +10.82 | +5.98 |
| ν=0.49 40×4 | −40.78 | −40.30 | +2.52 | +1.44 |

The overshoot is exact element behaviour, not a solver artefact. In the
one-element-thick bending mode the centroid dilatation is zero and ε_yy
cannot follow the linear Poisson contraction, so only the deviatoric
modulus resists: `(4/3)G` instead of `E'`. The predicted ratio
E'/((4/3)G) = 2.14 (ν=0.3) and 2.94 (ν=0.49) matched the removed
`-both` 10×1 (+113 %, +192 %), where the reduced shear adds no
bending resistance either.

For the plane-stress cantilever (the C80 study case, ν=0.3; tip error
20×2), Gauss gives −10.71 % (= Pro), `-yes` −2.43 %, Pro f=2 −2.75 %
and Pro f=10 +31 %. `-volumetric` is dropped on membrane groups.

Hourglass (the C80 study cases; Gauss = Pro gauss exactly):

| | Gauss | `-yes` | `-volumetric` | Pro f=2 | Pro f=10 |
|-|-------|--------|---------------|---------|----------|
| single quad4, hourglass load, max \|u\| | 3.47e-5 | 4.46e-5 | 6.69e-5 | 4.62e-5 | 1.16e-3 |
| 4×4 block, point load, mid-row zig-zag | 0.104 | 0.090 | 0.076 | 0.088 | 0.799 |

Isolated-element eigenvalues (independent reconstruction): quad4: 3
zero eigenvalues (the rigid modes) for Gauss, `-yes` and `-volumetric`.
hex8: 6 for Gauss and `-volumetric`, and 9 (6 rigid + 3 twist) for
`-yes`. `-volumetric` adds no spurious mode.

Tests (`validation-suite/test-2026`, registered in build_safe.sh):

| test | tolerance |
|------|-----------|
| `sri_volumetric_lame049_1` | u_r(a) 1e-5 = 0.5 %, σ_θ(a) 0.1 = 6 % |
| `sri_volumetric_lame04999_1` | as above |
| `sri_volumetric_hex8_lame1` | as above |

Each tolerance is set above the measured error against the analytic
value: u_r 2.6× (−0.19 %), σ_θ 1.9× (−3.1 %, from the nodal recovery at
the loaded edge). Gauss fails every target by an order of magnitude.
RED on the previous binary: rc 1, "I do not know : -volumetric". (The
C80.1 test `sri_both_cantilever1` was deleted with `-both` in C80.2.) `-yes` is covered by the
unchanged families `qsri` and `qsri3d`.

### Limitations / pending

- Linear elasticity only (the existing scope gate). B-bar for
  plasticity (plastic incompressibility) needs the reduced point to have
  its own volumetric state; not implemented, and an owner decision.
- Axisymmetric and plane-stress groups are refused (warning, volumetric
  bit dropped).
- The hex8 uses the 1-point centroid b in both the stiffness and the
  stored stress. For a distorted (non-parallelepiped) hex8, the Gauss
  sum Σ b_ip w detJ in the section post-processing differs slightly
  from 8 detJ_c b_c (calcul_force only).
- κ = m^T D m / 9 is exact for isotropic tangents. An anisotropic
  tangent still gets an algebraically exact split D_dev + κ m m^T, but
  D_dev is then not purely deviatoric.

## C80.2 — `-lambda` (Hughes 1977 λ/μ split) and removal of `-both` (2026-10-16)

Owner decision 2026-10-16: add `-lambda` (only the λ·m·mᵀ term of the
isotropic elastic matrix at the centroid point, the 2μ part at the full
Gauss rule; linear elastic materials), remove `-both` (nobody used it,
never better than `-volumetric`), keep `-yes` and `-volumetric`
unchanged. Gate (owner): `-lambda` measured on the C80.1 problems
before shipping; it may not be materially worse than `-volumetric` or
Gauss, nor show a spurious mode. The gate passed (tables below).

### Formulation

```
D = 2 mu I~ + lambda m m^T          (isotropic; I~ = identity with 1/2 on the shear entries)
K = sum_GaussIP B^T (D - lambda m m^T) B w detJ  +  lambda * volfac * W * detJ_c * b b^T
```

It is the C80.1 machinery with the coefficient κ replaced by λ:
`sigvec[0,3,5] -= λ·θ_ip` (D-c feedback), `new_sig[aa] += λ·(θ_c −
θ_ip)` (stored stress 2μ·ε_ip + λ·θ_c·I), `ddsdde_total[0,3,5 × 0,3,5]
-= λ`, and the reduced stiffness λ·W·detJ_c·volfac·b bᵀ joins
`stiffness_shear`. The energy is 2μ ε:ε at the Gauss points plus
λ·θ_c² at the centroid. The 2μ ε:ε part still holds the volumetric
share 2μ/3·(tr ε)² at the Gauss points; it does not lock because μ is
bounded as ν → 0.5 (measured: Lamé u_r(a) −0.24 % at ν = 0.49 and at
ν = 0.4999; the cantilever error does not grow from ν = 0.49 to 0.4999).
Positive definiteness: for κ = λ + 2μ/3 > 0 the element energy is
≥ (2μ/3 + λ)·θ_c² plus the deviatoric energy (Jensen on the IP
dilatations), so no extra zero mode appears; the eigenvalue check
confirms it.

Why it bends better than `-volumetric` on coarse meshes: in the
one-element-thick bending mode θ_c = 0, so `-volumetric` resists with
the deviatoric modulus only, while `-lambda` keeps 2μ/3·θ_ip² at the
Gauss points, which carries most of the missing Poisson coupling.

**Non-isotropic materials: refused.** λ is defined only for an
isotropic tangent. `materi()` takes λ as the mean of the 6 off-diagonal
entries of the normal 3×3 block of `ddsdde_total` and checks the block
is isotropic (equal diagonal entries, equal off-diagonal entries,
relative tolerance 1e-6 of the largest entry). Otherwise it stops with
"group_element_selective_reduced_integration -lambda needs an isotropic
elastic material ...; use -volumetric for this group. Element : n".
Defining an "effective λ" from an anisotropic D was rejected: the split
D − λ m mᵀ is no longer the λ/μ split, its accuracy is unmeasured, and
`-volumetric` already handles any tangent with κ = mᵀDm/9.

**Plane stress: dropped**, consistent with `-volumetric`. In a membrane
group ε_zz is condensed out, the effective λ* = 2λμ/(λ+2μ) stays
bounded (≤ 2μ) and there is no volumetric locking to remove;
`sri_active()` drops the `SRI_LAMBDA` bit, the group runs with the full
Gauss rule and the one-time warning is printed (measured: the
plane-stress cantilever equals Gauss to all digits). Plane-stress
bending is the job of `-yes`.

### Files and functions

- `tochnog.h` / `tochnog-mod.h`: the enum `BOTH` renamed `LAMBDA` IN
  PLACE (same position before `LAST_DUMMY`, so no enum value moved);
  `#define SRI_LAMBDA 4` next to `SRI_SHEAR` / `SRI_VOLUMETRIC`.
- `database.cc`: `name[LAMBDA] = "lambda"` (the name `both` is gone).
- `input.cc` (data-part reader, the "I do not know" branch): a
  dedicated error for `-both` on
  `group_element_selective_reduced_integration`, naming `-lambda` and
  `-volumetric`.
- `miscel.cc` `sri_active()`: accepts -no/-yes/-volumetric/-lambda
  (anything else `db_error`); `-lambda` returns `SRI_LAMBDA` unless the
  group is a membrane.
- `materi.cc` `materi()`: the requested mask and warning know
  `-lambda`; the C80.1 volumetric block and reduced-stiffness block run
  for `SRI_VOLUMETRIC | SRI_LAMBDA`; for `SRI_LAMBDA` the coefficient
  `sri_kappa` is λ from the isotropy-checked tangent.
- `calcul_force.cc`: unchanged. The section post-processing replaces
  the shear only for `SRI_SHEAR`; the stored λ·θ_c part is constant per
  element, so its Gauss-rule internal force equals the reduced one, as
  for `-volumetric`.

### Verification (measured 2026-10-16)

Independent checks: a pure-displacement Python solve of the same
element (2×2 Gauss D − λ m mᵀ + 1-point λ b bᵀ) reproduces the GNU
staggered results to 1e-10 relative on a distorted 20×2 cantilever
(ν = 0.3 / 0.49 / 0.4999) and to 1e-12 on a randomly distorted 3×2×2
hex8 block (ν = 0.3 / 0.4999). Gauss, `-yes` and `-volumetric` are
unchanged (the omen corpus A/B below).

Lamé cylinder (u_r(a) error % / σ_θ(a) node 1 error % / max deviation
of the IP mean stress % of p a²/(b²−a²)):

| ν, mesh | Gauss | `-yes` | `-volumetric` | `-lambda` |
|---------|-------|--------|---------------|-----------|
| 0.3, 8×16 | −0.33 / +5.7 / 31 | −0.26 / +6.4 / 32 | −0.17 / −3.1 / 4 | −0.22 / −0.2 / 13 |
| 0.49, 8×16 | −4.01 / +198 / 606 | −3.94 / +207 / 607 | −0.19 / −3.1 / 4 | −0.24 / −0.1 / 13 |
| 0.4999, 4×8 | −94.05 / +2145 / 6625 | −94.05 / +2322 / 6626 | −0.74 / −5.8 / 8 | −0.95 / −0.6 / 22 |
| 0.4999, 8×16 | −80.16 / +4196 / 12531 | −80.16 / +4370 / 12533 | −0.19 / −3.1 / 4 | −0.24 / −0.1 / 13 |
| 0.4999, 16×32 | −50.39 / +5697 / 16729 | −50.38 / +5815 / 16731 | −0.05 / −1.6 / 2 | −0.06 / −0.0 / 7 |

Pro `-reduced` (C80 study, 8×16): f=2 −0.25 / −3.04 / −75.18 % u_r at
ν = 0.3 / 0.49 / 0.4999, f=10 −0.01 / −0.13 / −10.80 %. The hex8 Lamé
(one layer, w = 0) equals the quad4 to all digits for every value.

Checkerboard: the `-lambda` IP mean-stress deviation (13 %) is the
intra-element variation μ·θ_ip of the 2μ share. It is ν-independent
(13 % at ν = 0.3 too, against 31 % for plain Gauss there) and O(h)
(22 → 13 → 7 %). The element-average pressure oscillation (|p_e − mean
of the edge neighbours|) is 0.1 % for both `-lambda` and `-volumetric`
against 109 % for Gauss at ν = 0.4999: no checkerboard.

Bending, tip deflection error % (plane strain; pure bending = constant
moment σ_xx = y, exact clamp; hex8 = one layer, w = 0, identical to the
quad4 to all digits):

| case | Gauss | `-yes` | `-volumetric` | `-lambda` |
|------|-------|--------|---------------|-----------|
| TG ν=0.3 10×1 / 20×2 / 40×4 / 80×8 | −36.28 / −12.42 / −3.44 / −0.88 | −18.32 / −5.28 / −1.38 / −0.35 | +22.08 / +4.70 / +1.14 / +0.28 | −4.86 / −1.25 / −0.32 / −0.08 |
| TG ν=0.49 | −91.23 / −70.61 / −40.78 / −15.36 | −91.09 / −70.24 / −40.30 / −15.09 | +67.15 / +10.82 / +2.52 / +0.62 | +30.30 / +5.98 / +1.44 / +0.36 |
| TG ν=0.4999 | −98.56 / −92.31 / −90.39 / −86.86 | same as Gauss | +70.42 / +11.14 / +2.59 / +0.64 | +32.85 / +6.37 / +1.53 / +0.38 |
| pure bending ν=0.3 | −36.47 / −12.51 / −3.46 / −0.89 | −18.35 / −5.30 / −1.38 / −0.35 | +22.42 / +4.78 / +1.16 / +0.29 | −4.76 / −1.23 / −0.31 / −0.08 |
| pure bending ν=0.49 | −92.24 / −72.67 / −41.91 / −15.70 | −92.09 / −72.29 / −41.42 / −15.43 | +67.90 / +11.04 / +2.56 / +0.63 | +30.65 / +6.11 / +1.47 / +0.36 |
| pure bending ν=0.4999 | −99.67 / −95.38 / −93.78 / −90.19 | same as Gauss | +71.22 / +11.38 / +2.63 / +0.65 | +33.22 / +6.51 / +1.56 / +0.39 |
| TG 20×2 distorted ±25 % ν=0.3 / 0.49 / 0.4999 | −19.33 / −61.61 / −90.12 | −14.78 / −61.22 / −90.12 | −8.40 / −0.72 / −0.30 | −12.25 / −3.79 / −3.32 |
| TG plane stress ν=0.3 10×1 / 20×2 / 40×4 | −32.46 / −10.71 / −2.92 | −9.09 / −2.43 / −0.62 | = Gauss (dropped) | = Gauss (dropped) |

Both selective volumetric values converge at O(h²) in bending
(factor ~4 per refinement); `-lambda` has 3-5× smaller errors at
ν = 0.3 and about half at ν ≥ 0.49 on regular meshes. On the distorted
mesh, `-lambda` is 3-4 points below `-volumetric`; this is not
locking (no growth from ν = 0.49 to 0.4999), it is the distortion
stiffening that `-volumetric` partly cancels with its regular-mesh
softness (+11 % on the regular 20×2).

Hourglass and eigenvalues:

| | Gauss | `-yes` | `-volumetric` | `-lambda` |
|-|-------|--------|---------------|-----------|
| single quad4, hourglass load, max \|u\| | 3.47e-5 | 4.46e-5 | 6.69e-5 | 5.20e-5 |
| 4×4 block, point load, mid-row zig-zag | 0.104 | 0.090 | 0.076 | 0.083 |
| zero eigenvalues, quad4 square / distorted | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 |
| zero eigenvalues, hex8 brick / distorted | 6 / 6 | 9 / 6 to 8 (near-singular, smallest 3e-3) | 6 / 6 | 6 / 6 |

No negative eigenvalue for any value; the smallest non-zero eigenvalue
of `-lambda` is above that of `-volumetric` (quad4 distorted ν=0.4999:
319 vs 275; hex8 distorted 25.9 vs 25.5).

Tests (`validation-suite/test-2026`, registered in build_safe.sh):

| test | target |
|------|--------|
| `sri_lambda_lame04999_1` | u_r(a) 1e-5 = 0.5 % (measured −0.24 %), σ_θ(a) 0.02 = 1.2 % (measured −0.14 %; `-volumetric` −3.1 % would fail) |
| `sri_lambda_cantilever1` | 10×1 ν=0.3 tip v(L, D/2) 6 % (measured −4.86 %; Gauss −36 %, `-volumetric` +22 %) |
| `sri_both_removed_error1` | expected error "integration -both was removed" |
| `sri_lambda_anisotropic_error1` | expected error "integration -lambda needs an isotropic elastic material" (transversely isotropic group) |

RED on the previous binary: the two `-lambda` targets rc 1 "I do not
know : -lambda"; the `-both` input ran with rc 0; the anisotropic input
stopped with "I do not know : -lambda" (not the expected message).

Suite and corpus (commit `5993499`): clean build; build_safe
TN_CHECK_ONLY 286 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0.
Omen corpus A/B against the C80.1 binary (977 runs): 0 node_dof
changes. Text-only differences: `bicg_csr_threads1` (threaded
round-off), `alternate1` (finished at 119.8 s against the 120 s timeout,
identical node_dof), `cylinder_buckle` (timeout in both); the removed
`sri_both_cantilever1` input now stops with the `-both` error, as
intended.

### Limitations / pending

- Linear isotropic elasticity only; plasticity stays with the planned
  `-volumetric` (B-bar) path.
- Plane stress and axisymmetry: dropped / refused, as `-volumetric`.
- The `expected_error_ok` helper of build_safe.sh passes the marker to
  `grep -qF` without `--`, so a marker that starts with `-` is read as
  a grep option; the markers of the C80.2 tests start with a word.

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
