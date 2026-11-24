# group_materi_elasti_shear_factor

## Description

`group_materi_elasti_shear_factor` (manual Professional 6.654) multiplies
the **shear stiffness** following from a specified young and poisson by a
`factor`:

> Specifying this record causes the shear stiffness following from a
> specified young and poisson to be multiplied with `factor`. This
> provides a convenient way to test in a calculation what the effect of
> low shear stresses is.

The normal-stiffness entries of the elasticity tensor are NOT touched:
only the shear blocks `C[i][j][k][l]` with `i != j` and `k != l` are
scaled (the Voigt diagonal shear entries (3,3), (4,4), (5,5) in 3D and
the (2,2) shear entry in 2D plane strain/plane stress, plus their
symmetric twins).

## Uso

```
group_type 0  -materi
group_materi_elasti_young 0  1000.0
group_materi_elasti_poisson 0  0.0
group_materi_elasti_shear_factor 0  2.0
```

A factor below 1 reduces the shear stiffness (e.g. 0.1 to model very low
shear stresses); `0.0` makes the element shear-soft completely.

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_elasti_shear_factor` | `factor` | Multiplication factor for the shear stiffness derived from young + poisson. |

## Validation

- `mshf.dat` / `mshf_nof.dat` (validation-suite/test-2014): single quad4
  in pure shear (bottom -velx, top +velx), 1 step (`eptxy = 0.03216`),
  `E = 1000`, `nu = 0` (`G = 500`).
  - factor `2.0`: `sigma_xy = 50`; no record: `sigma_xy = 33.33`;
    factor `0.0`: `sigma_xy = 0` EXACT (all the shear stress flows
    through the scaled entries — the discrimination is hard).
  - NOTE: the measured ratio is 1.5, not 2, because the post_point
    stress dof is a solved unknown of the coupled stress formulation
    (see the developer manual); the strain dof also shifts slightly
    (`eptxy 0.0322 -> 0.025` with factor 2). The direction of the effect
    (stiffer shear -> higher shear stress) and the exact zero at
    factor 0 are the validated physics.

## Notas

- Requires `materi_stress` and `materi_velocity` in the initialization
  part.
- Applies to the C and Cmem tensors after the young / young_polynomial /
  young_power blocks, so it combines with all constant and pressure-
  dependent young laws (but not with the volumetric-young / camclay /
  lade elastic laws, which are separate models).
