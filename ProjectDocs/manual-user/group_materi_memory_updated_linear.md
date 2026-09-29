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

- Rigid body rotations are neglected: the extra stresses are added to the
  previous stresses and the result is NOT rotated (manual Professional
  2.2.1). Until 2026-09-29 the GNU rotated the stress by the incremental
  rotation, which gave wrong stresses in shear (simple shear gamma = 1 with
  Drucker-Prager c = 0.1: sigxy 0.0693 instead of 0.11547, spurious
  sigxx/sigyy +-0.0924). Now identical to the Professional.
- `group_materi_elasti_compressibility` adds `tr(incremental elastic
  strain)/compressibility` per step. Until 2026-09-29 the GNU added the
  TOTAL strain every step unless `materi_strain_elasti` was initialized
  (two bars in series: -0.01 instead of the analytic -0.005).
- The rotation handling equals `-updated_without_rotation`; the strain
  measure (linear engineering strains) and the compressibility block equal
  `-total_linear`.
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
- Difference with the Professional on a moving mesh: the Professional's
  default mesh is `-fixed_in_space` (its manual 6.806), the GNU's default is
  `-follow_material` when `materi_velocity` is used without
  `materi_displacement`. Without a `mesh` record both give the engineering
  strain. With an explicit `mesh -follow_material` the Professional builds
  the increment on the current geometry (1D bar stretched 50%: strain
  0.40714, logarithmic-like) while the GNU still gives 0.5. Give the `mesh`
  record explicitly if the result must match in that case.
