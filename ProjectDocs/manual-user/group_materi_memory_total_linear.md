# group_materi_memory -total_linear

## Description

`group_materi_memory -total_linear` (manual Professional 6.685, theory
2.2.1) is the total Lagrange formulation with LINEAR engineering strains.
Deformations are measured from time 0 with the total deformation gradient
`F`; rigid-body rotations are neglected and the strain is the linear
engineering strain `0.5(F + F^T) - I`.

| Memory | Kinematics |
|--------|-----------|
| `-updated` | updated Lagrange, rate-of-deformation strains |
| `-updated_linear` | updated Lagrange, LINEAR engineering strains |
| `-total_linear` | total Lagrange, linear engineering strains |
| `-total` | total Lagrange, finite strains |

## Usage

```
group_type            0  -materi
group_materi_memory   0  -total_linear
```

## Behaviour

The stress at a new time point is the PREVIOUS stress plus the extra
stress caused by the INCREMENTAL strain (the difference of the linear
engineering strains between the two time points). This is important when
the initial stress state is prescribed with `control_reset_dof`: the
reset value is preserved by the next solve (the incremental strain is
zero when the displacements do not change).

## Notes

- This is the formulation used by the delete/excavation family
  (`control_mesh_delete_geometry` with a reset pre-stress): `delete3`
  (target `sigyy = 0.6`) and `mpc5` (target `sigyy = 0.5`) set an initial
  vertical stress with `control_reset_dof -sigyy` and then delete part of
  the mesh.
- The reset value survives only with the incremental update. Before this
  fix the GNU recomputed the stress from the total strain, discarding the
  reset (the delete tests gave `sigyy = 0`).
- The reported stress in these problems converges with the number of
  `control_timestep_iterations`; the Professional's default is 3, the GNU
  default is 1.
- `group_materi_elasti_compressibility` adds `tr(incremental elastic
  strain)/compressibility` per step (since 2026-09-29). Before, the GNU
  added the TOTAL strain every step unless `materi_strain_elasti` was
  initialized, which made the result depend on the time step. The
  Professional 02-08-2026 still shows that behaviour with `-total_linear`
  (two bars in series, analytic -0.005: Pro -0.0075 / -0.0099 / -0.009999
  for dt 1e-2 / 1e-3 / 1e-4; GNU -0.005 for every dt). This is a deliberate
  difference with the Professional.

## Differences with the Professional version

- **`-total_linear` + `group_materi_elasti_compressibility`** (deliberate
  divergence, project decision 2026-09-29; BACKLOG C9; SEGUIMIENTO
  Descarte). GNU: Volumetric increment from `tr(inc_epe)`: sigxx -0.005 for
  every time step (analytic, manual 2.2.1). Professional 02-08-2026:
  Time-step dependent: -0.0075 / -0.0099 / -0.009999 for dt 1e-2 / 1e-3 /
  1e-4 (keeps the legacy `tr(new_epe-old_epe)` increment). Status:
  Deliberate divergence, project owner (merge `93966ef`); GNU follows manual
  + analytic. Also Pro-suspect. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
