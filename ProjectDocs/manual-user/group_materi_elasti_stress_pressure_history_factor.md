# group_materi_elasti_stress_pressure_history_factor

## Description

`group_materi_elasti_stress_pressure_history_factor` (manual Professional
6.655) models a **different soil stiffness on first loading versus
unloading/reloading**: while the current pressure is SMALLER than the
largest pressure in history, the material is unloading or reloading and
the elastic stiffness is **multiplied with `factor`**; when the current
pressure is the new maximum, it becomes the maximum history pressure and
the stiffness is **not** multiplied. Typically `factor ~ 3`.

Requires the initia [`materi_stress_pressure_history`]
(materi_stress_pressure_history.md) (manual 4.50), which stores the
maximum of the absolute value of the pressure over time in the node_dof
records.

It can be combined with the young specified by
`group_materi_elasti_young` or the young calculated from
`group_materi_elasti_young_power` (and with
`group_materi_elasti_poisson_power` / `group_materi_elasti_shear_factor`).

## Uso

```
group_type 0  -materi
group_materi_elasti_young 0  1000.0
group_materi_elasti_poisson 0  0.3
group_materi_elasti_stress_pressure_history_factor 0  3.0
```

with, in the initialization part:

```
materi_velocity
materi_velocity_integrated
materi_strain_total
materi_stress
materi_stress_pressure_history
```

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_elasti_stress_pressure_history_factor` | `factor` | Stiffness multiplier applied while unloading/reloading (`|p| < max |p| history`). Typically ~3. |

## Física

The pressure is `p = -sig_mean` (positive in compression, same
convention as `group_materi_elasti_young_power`). During virgin
loading the pressure grows and each step becomes the new history
maximum -> the factor is NOT applied. During unloading or reloading the
pressure stays below the history maximum -> the factor IS applied: the
unloading/reloading tangent is stiffer, so the material recovers less
strain per stress drop (and shows the classic over-recovery to tension
when fully unloaded with a large factor).

## Validation

- `msph.dat`: oedometer with `E = 1000`, `nu = 0.3`, factor 3, two
  phases — loading (`vely = -0.002` for 0.2 s, peak
  `sigma_zz = -0.5385`, max `|p| = 0.3333` stored in the `sph` dof) and
  unloading (`vely = +0.002` for 0.1 s). With factor 3 the unloading
  tangent is `3*E` and the stress overshoots to
  `sigma_zz = +0.2692`; the `sph` dof stays at the peak `0.3333`.
- `msph_flat.dat`: same input with factor 1 (no effect): the material
  recovers elastically to `sigma_zz = -0.2692` (A/B discriminates the
  unloading stiffness); `sph` still `0.3333`.

## Notas

- The `sph` dof can be plotted / targeted with the basename `sph`
  (e.g. `target_item 1 -post_point_dof 1 -sph`).
- Requires `materi_stress`, `materi_velocity` and
  `materi_stress_pressure_history` in the initialization part.
