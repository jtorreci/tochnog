# materi_strain_total_kappa

## Description

`materi_strain_total_kappa` (manual Professional 4.45) adds a scalar dof (label `kapto`) to the
`node_dof` records: the maximum size of the total strain, `max sqrt(0.5 eps:eps)` as occurred in the history of the material
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
materi_strain_total_kappa
```

## Example

`validation-suite/test-2026/materi_strain_total_kappa1.dat`: one quad4
with the homogeneous strain eps_xx 0.004, eps_yy -0.01, unloaded to half
in the second step. The dofs keep the peak values: kappa 0.00761577,
compression -0.01, shear 0.00721110, tension 0.004 (analytic).

## Differences with the Professional version

The Professional 02-08-2026 adds the same dof (same label). On the
analytic test it gives the same kappa (history maximum of
`sqrt(0.5 a:a)`) and the same negative compression value at the peak,
but after unloading it reports the current compression and shear values
instead of the history value, and its tension kappa stays 0 (it copies the
shear value when both are initialised). GNU follows the manual ("as
occurred in history") and takes the negative compression sign from the
Professional (DIVERGENCES STRAIN-KAPPA-HISTORY, BACKLOG C108).
