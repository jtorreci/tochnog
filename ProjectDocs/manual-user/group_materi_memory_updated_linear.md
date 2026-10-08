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

`-updated_linear` is also the DEFAULT memory when a group has no
`group_materi_memory` record (since 2026-10-04, as the Professional); see
[group_materi_memory_default](group_materi_memory_default.md); with
`materi_displacement` the default is `-total_linear` instead (F3).

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
- The strain increment is the linear engineering strain of the velocity
  gradient on the CURRENT geometry (manual Professional 2.2.1). With the
  default `-fixed_in_space` mesh the geometry does not move, so the
  accumulated total is the engineering strain (`eptyy = -0.3` for a 30%
  compression, hypo7/8/9). With an explicit `mesh -follow_material` the
  increments are summed on the moving geometry (bar stretched by 50%:
  0.407136, as the Professional).
- Until F3 (2026-10-07) the GNU built the increment of an explicit
  `-updated_linear` on the reference geometry (`node_start_refined`) to
  emulate a fixed mesh under the old GNU `-follow_material` default; that
  path was removed when the default mesh became `-fixed_in_space`.
- Self-checking inputs: `validation-suite/test-2026/materi_memory_updated_linear_current1.dat`
  (bar, sigxx 0.407136) and `materi_memory_updated_linear_current2.dat`
  (column, 100 steps, node 3 velix 0.52093 within 1.e-3).

## Differences with the Professional version

- **explicit `-updated_linear` on a `mesh -follow_material`** (C19):
  RESOLVED in F3 (2026-10-07). The increment is built on the current
  geometry, as the Professional (bar +50%: 0.407136 in both).
- **`-updated_linear` together with `materi_displacement`** (known
  difference, documented in F3 2026-10-07). GNU: accepted (the
  coordinates are updated, elem.cc). Professional 02-08-2026: rejects the
  combination. Status: Open, LOW; GNU stays permissive. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **explicit `-follow_material`: in-step geometry** (known difference).
  GNU builds the increment on the END-of-step geometry, the Professional
  on the START-of-step geometry: 2-quad4 column sheared by 1, node 3, one
  step GNU 0.47222 / Pro 0.5, 100 steps GNU 0.52048 / Pro 0.52093.
  Status: Open, LOW. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
