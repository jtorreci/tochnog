# group_materi_elasti_c / group_materi_elasti_c_direction

## Description

`group_materi_elasti_c` (manual Professional 6.642) gives the 81 values of
the linear material stiffness C_ijkl directly, with
sigma_ij = C_ijkl eps_kl (tensor strains). The sequence of stresses and
strains is xx, xy, xz, yx, yy, yz, zx, zy, zz; the values are given row by
row, 9 values per row. Useful for any anisotropic linear material
(orthotropic rock, layered soils, homogenized composites).

The stiffness is added to the other elastic laws of the group (for example
`group_materi_elasti_young` / `_poisson`), as the Professional does.

`group_materi_elasti_c_direction` (6.643) gives the local axes dir_0,
dir_1, dir_2 (global components, 9 values) in which C is specified. Each
axis is normalized. Default: the global axes.

## Syntax

```
group_materi_elasti_c index c_xxxx c_xxxy ... c_zzzz      (81 values)
group_materi_elasti_c_direction index d0x d0y d0z d1x d1y d1z d2x d2y d2z
```

For an isotropic material C_ijkl = lambda d_ij d_kl + G (d_ik d_jl + d_il d_jk),
so the shear entries C_xyxy = C_xyyx = G.

## Example

```
group_materi_elasti_c 0
  1346.15 0 0 0 576.923 0 0 0 576.923
  0 384.615 0 384.615 0 0 0 0 0
  ...
group_materi_elasti_c_direction 0  0 1 0  -1 0 0  0 0 1
```

Tests: `validation-suite/test-2026/group_materi_elasti_c1` (isotropic C =
young / poisson), `_c2` (shear), `_c3` (added to young / poisson),
`group_materi_elasti_c_direction1`-`3`, `_c_error1` (membrane),
`_c_direction_error1` (zero axis).

## Differences with the Professional version

- A zero axis in `group_materi_elasti_c_direction` is an error; the
  Professional 02-08-2026 runs with a degenerate stiffness
  (DIVERGENCES ELASTI-C-ZERO-DIRECTION).
- Not available with `group_materi_membrane -yes` (error), as the
  Professional.
