# group_materi_memory -updated_linear

## Description

`group_materi_memory -updated_linear` (manual Professional 6.685) is
the updated Lagrange formulation with LINEAR strains (small deformation
theory on the current configuration). The GNU treated it as an unknown
memory type (the corpus hypoplasticity tests, e.g. hypo7/8, use it).

| Memory | Kinematics |
|--------|-----------|
| `-updated` | updated Lagrange, rate-of-deformation strains |
| `-updated_linear` | updated Lagrange, LINEAR engineering strains (small deformations) |
| `-total_linear` | total Lagrange, linear strains |
| `-total` | total Lagrange, finite strains |

## Usage

```
group_type            0  -materi
group_materi_memory   0  -updated_linear
```

## Notes

- In the GNU implementation `-updated_linear` behaves like
  `-updated_without_rotation` for the rotation handling (identity
  rotation) and like `-total_linear` for the strain measure (linear
  engineering strains) and the compressibility block.
- hypo7/hypo8 of the corpus use it (run and converge close to the
  Professional targets).
