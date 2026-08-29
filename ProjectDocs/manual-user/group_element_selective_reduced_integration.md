# group_element_selective_reduced_integration

## Description

Selective reduced integration (SRI, see Hughes, "The Finite Element
Method", chapter on reduced integration) for the bilinear quad4 (2D)
and the trilinear hex8 (3D), opt-in. It fixes the CLASSIC shear
locking / parasitic shear of the bilinear elements in bending:

- The bilinear quad4 / trilinear hex8 cannot represent the quadratic
  transverse displacement `u_y ~ -κ·x²/2` of a beam in bending. When
  forced to bend they generate a spurious shear strain `γ_xy ≠ 0`
  (parasitic shear), which makes the bending response TOO STIFF: the
  element under-estimates moments and displacements (measured in this
  GNU: a cantilever with one quad4 in the thickness gives 0.23× of the
  analytic section moment; the 3D hex8 cantilever gives 0.221× of the
  Euler-Bernoulli deflection).
- With SRI the shear terms of the constitutive matrix (γ_xy for the
  quad4; γ_xy, γ_xz, γ_yz for the hex8) are integrated with 1 Gauss
  point at the element centroid while the normal terms keep the full
  Gauss rule (2×2 / 2×2×2). The parasitic shear vanishes at the
  centroid for the bending mode, so the spurious shear energy
  disappears from the stiffness (measured 3D cantilever: the
  deflection improves from 0.221× to 0.897× of the Euler-Bernoulli
  value).

The keyword is an opt-in: without it the element keeps the default
behavior exactly.

## Usage

```
group_element_selective_reduced_integration <element_group> -yes
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `element_group` | Element group, see `element_group`. |

## Scope

- 2D bilinear quad4 and 3D trilinear hex8.
- LINEAR ELASTICITY only: with plasticity, damage, maxwell, large
  displacement (`materi_displacement`) or axisymmetric groups the
  keyword is ignored with a one-time warning.
- Default OFF: groups without the keyword are byte-identical to the
  previous behavior.
- KNOWN LIMITATION (hex8, measured 2026-08-29): the shear-only SRI
  hex8 retains zero-energy modes (3 twist modes of the isolated
  element; section-warping modes of a mesh). The load-orthogonal
  cantilever solves correctly (0.897× of the Euler-Bernoulli
  deflection vs 0.221× locked), but configurations whose boundary
  conditions do NOT suppress the warping modes (free lateral faces,
  simply-supported beams) are numerically at risk. The 2D quad4 SRI
  does not have this limitation. See the developer manual for the
  eigenvalue evidence.

## Example

```
group_element_selective_reduced_integration 0  -yes
```

With a 2D cantilever modelled with one quad4 in the thickness (plane
stress), the fix C/D of the staggered scheme recovers the classic
Hughes result: the section moment at the clamp = 0.9375·P·L (measured,
family `qsri`, 2026-08-28). With a 3D cantilever modelled with one
hex8 in the section (family `qsri3d`, 2026-08-29), the tip deflection
improves from 0.221× to 0.897× of the Euler-Bernoulli value and the
axial bending stress at the clamp from 0.26× to 1.06× of the analytic
value.
