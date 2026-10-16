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
  (`-volumetric`, added 2026-10-16). When the Poisson ratio is close to
  0.5, the full rule enforces the incompressibility at every
  integration point. The element cannot satisfy this, so it becomes far
  too stiff, and the pressure shows a checkerboard. Measured on a thick
  cylinder in plane strain at ν = 0.4999, the full Gauss rule gives the
  inner radial displacement 80 % too small. With `-volumetric`, the
  volumetric part of the material matrix (bulk modulus times the
  dilatation) uses 1 point at the centroid. This is selective
  volumetric integration, also called B-bar. The deviatoric part keeps
  the full Gauss rule. The error drops to 0.19 %, and the pressure is
  constant per element, so there is no checkerboard.

`-both` combines the two.

## Usage

```
group_element_selective_reduced_integration <element_group> <value>
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `element_group` | Element group, see `element_group`. |
| `value` | `-no` (default behaviour) / `-yes` (shear terms at 1 point) / `-volumetric` (volumetric part at 1 point) / `-both` (shear terms and volumetric part at 1 point). Any other value is an input error. |

In every case, the rest of the material matrix uses the full Gauss rule
(2×2 / 2×2×2). The keyword replaces the group's
`group_integration_method` for these elements, including `-reduced`.

## Scope

- 2D bilinear quad4 and 3D trilinear hex8.
- LINEAR ELASTICITY only. With plasticity, damage, maxwell, large
  displacement (`materi_displacement`) or axisymmetric groups, the
  keyword is ignored and a one-time warning is printed.
- The volumetric part (`-volumetric`, `-both`) needs plane strain or
  3D. A plane-stress group (`group_materi_membrane -yes`) has no
  volumetric locking, because ε_zz is free. On such a group the
  volumetric part is dropped with a one-time warning: `-volumetric`
  falls back to the default rule, and `-both` acts as `-yes`.
- Default OFF: groups without the keyword, or with `-yes`, give
  byte-identical results to the previous behaviour.

## Which value to use

| Situation | Value |
|-----------|-------|
| Thin structure in bending, ν well below 0.5, 1-2 elements through the thickness | `-yes` |
| Near-incompressible material (ν ≥ 0.45, undrained soil, rubber-like), plane strain or 3D | `-volumetric` |
| Near-incompressible material in bending | `-volumetric` first. Use `-both` only with at least 4 elements through the thickness (see below). |

Measured accuracy against analytic solutions (plane strain, quad4,
error of the displacement):

| Case | Gauss (full) | `-yes` | `-volumetric` | `-both` | Pro `-reduced` f=2 | Pro `-reduced` f=10 |
|------|--------------|--------|---------------|---------|--------------------|---------------------|
| Lamé cylinder 8×16, ν=0.3 | −0.33 % | −0.26 % | −0.17 % | −0.10 % | −0.25 % | −0.01 % |
| Lamé cylinder 8×16, ν=0.49 | −4.0 % | −3.9 % | **−0.19 %** | −0.11 % | −3.0 % | −0.13 % |
| Lamé cylinder 8×16, ν=0.4999 | −80 % | −80 % | **−0.19 %** | −0.11 % | −75 % | −11 % |
| Lamé cylinder 4×8, ν=0.4999 | −94 % | −94 % | −0.74 % | −0.42 % | −92 % | −32 % |
| pressure checkerboard (max deviation of the IP mean stress, 8×16, ν=0.4999) | 12531 % | 12533 % | **4 %** | 4 % | 13493 % | 9401 % |
| Cantilever 10×1, ν=0.3 | −36 % | −18 % | +22 % | +113 % | — | — |
| Cantilever 20×2, ν=0.3 | −12 % | −5.3 % | +4.7 % | +15 % | — | — |
| Cantilever 40×4, ν=0.49 | −41 % | −40 % | +2.5 % | +4.2 % | — | — |
| Cantilever 20×2, plane stress, ν=0.3 | −10.7 % | **−2.4 %** | (= Gauss, dropped) | (= `-yes`) | −2.8 % | +31 % |

The hex8 Lamé cylinder (one layer, w = 0) gives the same numbers as the
quad4 to all digits.

Notes:

- Unlike `-reduced` with a large factor, the selective rules have **no
  hourglass modes** on the quad4. A single quad4 under a
  self-equilibrated hourglass load deforms 1.9× (`-volumetric`) and
  3.4× (`-both`) as much as with Gauss, against 33× for Pro `-reduced`
  f=10. Under a point load on a 4×4 block, the zig-zag of the middle
  row is 0.076 (`-volumetric`) and 0.085 (`-both`), against 0.104 for
  Gauss and 0.80 for Pro f=10.
- `-volumetric` and `-both` are SOFT in bending on meshes that are
  coarse through the thickness. With one element in the thickness, the
  bending mode has zero dilatation at the centroid, and the linear
  Poisson contraction across the thickness cannot be represented. So
  only the deviatoric part resists the bending strain. The cantilever
  converges from above (+113 % → +15 % → +3.4 % for `-both` at
  10×1 → 20×2 → 40×4, ν=0.3). Adding the reduced shear (`-both`) makes
  it softer still. For a bending-dominated group with ν well below 0.5,
  `-yes` is the better choice.
- KNOWN LIMITATION (hex8, measured 2026-08-29): the shear-only SRI hex8
  (`-yes`, and therefore also `-both`) keeps zero-energy modes: 3 twist
  modes of the isolated element and section-warping modes of a mesh.
  Configurations whose boundary conditions do not suppress the warping
  modes (free lateral faces, simply-supported beams) are numerically at
  risk. `-volumetric` on the hex8 has no zero-energy mode besides the 6
  rigid modes. The 2D quad4 has none for any value.

## Example

```
group_element_selective_reduced_integration 0  -volumetric
```

A quarter of a thick cylinder (a=1, b=2, internal pressure p=1,
E=1000, ν=0.4999, plane strain, 8×16 quad4) gives the inner radial
displacement within 0.19 % of the Lamé solution (Gauss: −80 %) and
σ_θ(a) = 1.615 against 5/3 (Gauss: 71.6). See the tests
`validation-suite/test-2026/sri_volumetric_lame04999_1.dat`,
`sri_volumetric_lame049_1.dat`, `sri_volumetric_hex8_lame1.dat` and
`sri_both_cantilever1.dat`.

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
ν = 0.5. The selective `-volumetric` removes the volumetric locking
without hourglass modes.
