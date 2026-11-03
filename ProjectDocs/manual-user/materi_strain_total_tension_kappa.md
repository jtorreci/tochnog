# materi_strain_total_tension_kappa

## Description

`materi_strain_total_tension_kappa` (manual Professional 4.48) adds a scalar dof (label `eptkt`) to the
`node_dof` records: the maximum principal tension total strain, `max(eps_max, 0)` as occurred in the history of the material
point. The value never decreases. When the material unloads, the dof keeps
its peak value. The manual uses these dofs for softening through a
`dependency_diagram`.

"Size" means the same measure as `materi_plasti_kappa`, `sqrt(0.5 a:a)`.
The strain is the total strain of `materi_strain_total`, which must be
initialised too (otherwise: "Error: materi_strain_total should be
initialized.").

## Usage

```
materi_strain_total
materi_strain_total_tension_kappa
```

## Example

`validation-suite/test-2026/materi_strain_total_kappa1.dat`: one quad4
with the homogeneous strain eps_xx 0.004, eps_yy -0.01, unloaded to half
in the second step. The dofs keep the peak values: kappa 0.00761577,
compression 0.01, shear 0.00721110, tension 0.004 (analytic).

## Differences with the Professional version

The Professional 02-08-2026 accepts this initialisation but adds no dof.
GNU follows the manual (owner decision 2026-11-03, DIVERGENCES
STRAIN-PLASTI-MODEL-DOFS). The size measure `sqrt(0.5 a:a)` and the
positive sign of the compression kappa are GNU choices (the manual does
not define them).
