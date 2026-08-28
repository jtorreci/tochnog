# group_element_selective_reduced_integration

## Description

Selective reduced integration (SRI, see Hughes, "The Finite Element
Method", chapter on reduced integration) for the bilinear quad4,
opt-in. It fixes the CLASSIC shear locking / parasitic shear of the
quad4 in bending:

- The bilinear quad4 cannot represent the quadratic transverse
  displacement `u_y ~ -κ·x²/2` of a beam in bending. When forced to
  bend it generates a spurious shear strain `γ_xy ≠ 0` (parasitic
  shear), which makes the bending response TOO STIFF: the element
  under-estimates moments and displacements (measured in this GNU:
  a cantilever with one quad4 in the thickness gives 0.23× of the
  analytic section moment).
- With SRI the shear term of the constitutive matrix is integrated
  with 1 Gauss point at the element centroid while the normal terms
  keep the full 2×2 Gauss rule. The parasitic shear vanishes at the
  centroid for the bending mode, so the spurious shear energy
  disappears from the stiffness, and the normal terms keep the
  element stable (no hourglass modes; the 3 rigid modes stay exact
  zero modes).

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

- 2D bilinear quad4 only (hex8 SRI in 3D is future work).
- LINEAR ELASTICITY only: with plasticity, damage, maxwell, large
  displacement (`materi_displacement`) or axisymmetric groups the
  keyword is ignored with a one-time warning.
- Default OFF: groups without the keyword are byte-identical to the
  previous behavior.

## Example

```
group_element_selective_reduced_integration 0  -yes
```

With a cantilever modelled with one quad4 in the thickness (plane
stress), the section moment at the clamp improves from 0.23× (locked)
to 0.31× of P·(8−x) and the section shear from 0.74× to 0.65× of P
(measured, family `qsri`). Note on the solver: the classic textbook
result (moment ≈ P·(8−x)) requires a reliable solve of the u-system;
the mixed u-σ Bi-CG of this GNU is numerically limited for the
1-element-in-thickness quad4 (see the developer manual; the exact
reference values are 99.3% deflection / 93.75% moment).
