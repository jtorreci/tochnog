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
- The accumulated strain is referred to the UNDEFORMED (reference) mesh,
  so the reported total is the LINEAR engineering strain
  (`eptyy = -0.3` for a 30% compression), matching the Professional.
  This holds even with the default velocity-integrated mesh, where the
  nodes follow the material: the strain increment is built with the
  reference-geometry shape derivatives (`node_start_refined`), while the
  reported mesh geometry still follows the material.
- Only `-updated_linear` uses this reference-geometry path; the other
  memory options (`-updated`, `-updated_without_rotation`,
  `-updated_area`, `-total`, `-total_linear`) are unchanged.
- hypo7/hypo8/hypo9 of the corpus use it (the strain measure now matches
  the Professional exactly; the final stress still differs by the kernel
  route, e.g. hypo7 sigyy -235.6 vs the Professional -231.8).
