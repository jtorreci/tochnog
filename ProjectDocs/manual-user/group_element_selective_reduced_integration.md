# group_element_selective_reduced_integration

**GNU-only extension.** The Professional manual and the Professional
binary have nothing equivalent. The closest Professional feature is
`group_integration_method -reduced`, which reduces the rule for the
whole stress (see the comparison below).

## Description

Selective reduced integration (SRI, see Hughes, "The Finite Element
Method", chapter on reduced and selective integration) for the bilinear
quad4 (2D) and the trilinear hex8 (3D). It is opt-in and never the
default. It removes two classic kinds of locking:

- **Shear locking** in bending (`-yes`). The bilinear quad4 and the
  trilinear hex8 cannot represent the quadratic transverse displacement
  `u_y ~ -κ·x²/2` of a beam in bending. When they are forced to bend,
  they generate a spurious shear strain `γ_xy ≠ 0` (parasitic shear),
  so the bending response is TOO STIFF and the element under-estimates
  moments and displacements. Measured in this GNU: a cantilever with one
  quad4 in the thickness gives 0.23× of the analytic section moment, and
  the 3D hex8 cantilever gives 0.221× of the Euler-Bernoulli deflection.
  With `-yes`, the shear terms of the constitutive matrix (γ_xy for the
  quad4; γ_xy, γ_xz and γ_yz for the hex8) use 1 Gauss point at the
  element centroid, and the normal terms keep the full Gauss rule
  (2×2 / 2×2×2).
- **Volumetric locking** of near-incompressible materials
  (`-volumetric` and `-lambda`, added 2026-10-16). When the Poisson
  ratio is close to 0.5, the full rule enforces the incompressibility
  at every integration point. The element cannot satisfy this, so it
  becomes far too stiff, and the pressure shows a checkerboard.
  Measured on a thick cylinder in plane strain at ν = 0.4999, the full
  Gauss rule gives the inner radial displacement 80 % too small.
  - With `-volumetric`, the volumetric part of the material matrix
    (bulk modulus κ times the dilatation) uses 1 point at the
    centroid. This is selective volumetric integration, also called
    B-bar. The deviatoric part keeps the full Gauss rule. The error
    drops to 0.19 %, and the pressure is constant per element.
  - With `-lambda` (Hughes 1977), only the Lamé-λ part of the isotropic
    elastic matrix, λ·(tr ε)², uses 1 point at the centroid. The 2μ
    part keeps the full Gauss rule, including its own volumetric
    share 2μ/3·(tr ε)². That share does not lock, because μ stays
    bounded when ν → 0.5. The error is 0.24 % at ν = 0.4999 (the same
    as at ν = 0.49), and the element-to-element pressure shows no
    checkerboard. In bending on coarse meshes `-lambda` is much more
    accurate than `-volumetric` (see the table).

The value `-both` (C80.1: shear terms and volumetric part at 1 point)
was removed on 2026-10-16 (C80.2). It was never more accurate than
`-volumetric`, and it was very soft in coarse bending. An input with
`-both` stops with the error "group_element_selective_reduced_integration
-both was removed; use -lambda (elastic, bending-dominated) or
-volumetric (general) instead."

## Usage

```
group_element_selective_reduced_integration <element_group> <value>
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `element_group` | Element group, see `element_group`. |
| `value` | `-no` (default behaviour) / `-yes` (shear terms at 1 point) / `-volumetric` (volumetric part κ·m·mᵀ at 1 point) / `-lambda` (the λ·m·mᵀ part at 1 point; isotropic elasticity only). Any other value is an input error. |

In every case, the rest of the material matrix uses the full Gauss rule
(2×2 / 2×2×2). The keyword replaces the group's
`group_integration_method` for these elements, including `-reduced`.

## Scope

- 2D bilinear quad4 and 3D trilinear hex8.
- LINEAR ELASTICITY only. With plasticity, damage, maxwell, large
  displacement (`materi_displacement`) or axisymmetric groups, the
  keyword is ignored and a one-time warning is printed.
- `-lambda` needs an ISOTROPIC elastic material, because λ is defined
  only for an isotropic material matrix. On another elastic material
  (for example `group_materi_elasti_transverse_isotropy`) the
  calculation stops with an error that names `-volumetric` as the
  alternative. `-volumetric` uses κ = mᵀ·D·m/9 and accepts any
  elastic material matrix.
- The volumetric part (`-volumetric`, `-lambda`) needs plane strain or
  3D. A plane-stress group (`group_materi_membrane -yes`) has no
  volumetric locking: ε_zz is free and the effective λ is
  2λμ/(λ+2μ), which stays bounded as ν → 0.5. On such a group the
  value is dropped with a one-time warning, and the group uses the
  default (full Gauss) rule. Use `-yes` for plane-stress bending.
- Default OFF: groups without the keyword, or with `-yes`, give
  byte-identical results to the previous behaviour.

## Which value to use

| Situation | Value |
|-----------|-------|
| Thin structure in bending, ν well below 0.5, 1-2 elements through the thickness, plane stress | `-yes` |
| Linear elastic isotropic material, bending-dominated or near-incompressible, plane strain or 3D | `-lambda` |
| Near-incompressible material that is not isotropic, or a general / mixed problem | `-volumetric` |
| Future plasticity (B-bar for elasto-plastic groups, not implemented yet) | `-volumetric` is the planned path |

Measured accuracy against analytic solutions (plane strain, quad4,
error of the displacement; the plane-stress row is a membrane group):

| Case | Gauss (full) | `-yes` | `-volumetric` | `-lambda` | Pro `-reduced` f=2 | Pro `-reduced` f=10 |
|------|--------------|--------|---------------|-----------|--------------------|---------------------|
| Lamé cylinder 8×16, ν=0.3 | −0.33 % | −0.26 % | −0.17 % | −0.22 % | −0.25 % | −0.01 % |
| Lamé cylinder 8×16, ν=0.49 | −4.0 % | −3.9 % | **−0.19 %** | **−0.24 %** | −3.0 % | −0.13 % |
| Lamé cylinder 8×16, ν=0.4999 | −80 % | −80 % | **−0.19 %** | **−0.24 %** | −75 % | −11 % |
| Lamé cylinder 4×8, ν=0.4999 | −94 % | −94 % | −0.74 % | −0.95 % | −92 % | −32 % |
| σ_θ(a) at the inner node, 8×16, ν=0.4999 | +4196 % | +4370 % | −3.1 % | **−0.14 %** | | |
| element pressure oscillation (8×16, ν=0.4999) | 109 % | 109 % | 0.1 % | 0.1 % | | |
| max deviation of the IP mean stress (8×16, ν=0.4999) | 12531 % | 12533 % | 4 % | 13 % | 13493 % | 9401 % |
| Cantilever 10×1, ν=0.3 | −36 % | −18 % | +22 % | **−4.9 %** | | |
| Cantilever 20×2, ν=0.3 | −12 % | −5.3 % | +4.7 % | **−1.3 %** | | |
| Cantilever 10×1, ν=0.49 | −91 % | −91 % | +67 % | +30 % | | |
| Cantilever 20×2, ν=0.49 | −71 % | −70 % | +11 % | +6.0 % | | |
| Cantilever 40×4, ν=0.49 | −41 % | −40 % | +2.5 % | **+1.4 %** | | |
| Cantilever 40×4, ν=0.4999 | −90 % | −90 % | +2.6 % | +1.5 % | | |
| Pure bending 10×1, ν=0.3 | −36 % | −18 % | +22 % | **−4.8 %** | | |
| Cantilever 20×2 distorted (±25 %), ν=0.3 | −19 % | −15 % | −8.4 % | −12 % | | |
| Cantilever 20×2 distorted (±25 %), ν=0.4999 | −90 % | −90 % | −0.3 % | −3.3 % | | |
| Cantilever 20×2, plane stress, ν=0.3 | −10.7 % | **−2.4 %** | (= Gauss, dropped) | (= Gauss, dropped) | −2.8 % | +31 % |

The hex8 versions (one layer of hex8 with w = 0: Lamé cylinder and
cantilever) give the same numbers as the quad4 to all digits.

Notes:

- The "max deviation of the IP mean stress" of `-lambda` (13 %) is the
  variation INSIDE each element of the 2μ share of the pressure. It
  does not grow with ν (13 % at ν = 0.3 as well, where plain Gauss
  already shows 31 %) and it halves with each mesh refinement
  (22 % → 13 % → 7 %). The element-to-element oscillation, which is
  what a checkerboard is, is 0.1 %, the same as `-volumetric`.
- On distorted meshes both `-volumetric` and `-lambda` lose some
  accuracy in bending. `-volumetric` looks better there only because
  its regular-mesh softness (+11 % on the regular 20×2) cancels the
  stiffening caused by the distortion. `-lambda` shows no locking
  trend: −3.8 % at ν = 0.49 and −3.3 % at ν = 0.4999.
- Convergence in mesh size (tip deflection, plane strain cantilever,
  10×1 → 20×2 → 40×4 → 80×8): `-lambda` −4.86 → −1.25 → −0.32 →
  −0.08 % at ν = 0.3, +30 → +6.0 → +1.4 → +0.36 % at ν = 0.49;
  `-volumetric` +22 → +4.7 → +1.1 → +0.28 % at ν = 0.3. Lamé u_r(a),
  4×8 → 8×16 → 16×32, ν = 0.4999: `-lambda` −0.95 → −0.24 → −0.06 %,
  `-volumetric` −0.74 → −0.19 → −0.05 %.
- Unlike `-reduced` with a large factor, the selective rules have **no
  hourglass modes** on the quad4. A single quad4 under a
  self-equilibrated hourglass load deforms 1.5× (`-lambda`) and 1.9×
  (`-volumetric`) as much as with Gauss, against 33× for Pro
  `-reduced` f=10. Under a point load on a 4×4 block, the zig-zag of
  the middle row is 0.083 (`-lambda`) and 0.076 (`-volumetric`),
  against 0.104 for Gauss and 0.80 for Pro f=10. The element
  stiffness of `-lambda` and `-volumetric` has only the rigid-body
  zero eigenvalues (3 for the quad4, 6 for the hex8), also on
  distorted elements.
- `-volumetric` is SOFT in bending on meshes that are coarse through
  the thickness. With one element in the thickness, the bending mode
  has zero dilatation at the centroid, and only the deviatoric part
  resists the bending strain. `-lambda` keeps the 2μ/3·(tr ε)² share
  at the full rule, which stiffens this mode back towards the exact
  value.
- KNOWN LIMITATION (hex8, measured 2026-08-29): the shear-only SRI hex8
  (`-yes`) keeps zero-energy modes: 3 twist modes of the isolated
  element and section-warping modes of a mesh. Configurations whose
  boundary conditions do not suppress the warping modes (free lateral
  faces, simply-supported beams) are numerically at risk.
  `-volumetric` and `-lambda` on the hex8 have no zero-energy mode
  besides the 6 rigid modes. The 2D quad4 has none for any value.

## Example

```
group_element_selective_reduced_integration 0  -lambda
```

A plane-strain cantilever (L=10, D=1, E=1000, ν=0.3, parabolic tip
shear P=1) with ONE row of 10 quad4 gives the tip deflection 4.9 % below
the Timoshenko-Goodier solution (Gauss: −36 %, `-volumetric`: +22 %).
A quarter of a thick cylinder (a=1, b=2, internal pressure p=1,
E=1000, ν=0.4999, plane strain, 8×16 quad4) gives the inner radial
displacement within 0.24 % of the Lamé solution (Gauss: −80 %) and
σ_θ(a) = 1.6644 against 5/3. See the tests
`validation-suite/test-2026/sri_lambda_cantilever1.dat` and
`sri_lambda_lame04999_1.dat`.

```
group_element_selective_reduced_integration 0  -volumetric
```

The same cylinder gives the inner radial displacement within 0.19 %
of the Lamé solution and σ_θ(a) = 1.615. See the tests
`sri_volumetric_lame04999_1.dat`, `sri_volumetric_lame049_1.dat` and
`sri_volumetric_hex8_lame1.dat`.

```
group_element_selective_reduced_integration 0  -yes
```

For a 2D cantilever modelled with one quad4 in the thickness (plane
stress), fix C/D of the staggered scheme recovers the classic Hughes
result: the section moment at the clamp is 0.9375·P·L (measured, family
`qsri`, 2026-08-28). For a 3D cantilever modelled with one hex8 in the
section (family `qsri3d`, 2026-08-29), the tip deflection improves from
0.221× to 0.897× of the Euler-Bernoulli value, and the axial bending
stress at the clamp improves from 0.26× to 1.06× of the analytic value.

## Differences with the Professional version

The Professional has no selective integration; this keyword is a GNU
extension. The Professional's `group_integration_method -reduced`
(implemented in the GNU as well, C80) reduces the rule for the WHOLE
stress with one factor. The table above shows the trade-off: f=10
removes most of the volumetric locking but overshoots bending and
allows hourglassing, while f=2 is good in bending but still locks near
ν = 0.5. The selective `-volumetric` and `-lambda` remove the
volumetric locking without hourglass modes.
