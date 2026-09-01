# materi_plasti_kappa_shear

## Description

`materi_plasti_kappa_shear` (manual Professional 4.25) is an
initialization option that adds a scalar dof `kapsh` to the node_dof
records containing the **size of the SHEAR part of the plastic strain**:

```
kappa_shear = int sqrt( 0.5 * dev(eps_p) : dev(eps_p) ) dt
```

i.e. the accumulated deviatoric plastic strain (the "equivalent plastic
shear strain"), independent of the total `kappa` dof
(`materi_plasti_kappa`, manual 4.24) which also includes the volumetric
part.

It is the hardening variable of the Mohr-Coulomb hardening-softening
family: `group_materi_plasti_mohr_coul_hardening_softening` and
`group_materi_plasti_mohr_coul_direct_hardening_softening` interpolate
the strength parameters linearly with `kappa_shear / kappa_shear_crit`
(manual Professional 6.731).

## Usage

In the initialization part:

```
materi_velocity
materi_displacement
materi_strain_total
materi_strain_plasti
materi_stress
materi_plasti_kappa
materi_plasti_kappa_shear
end_initia
```

The dof can be reset (e.g. when a construction layer is activated) and
used in targets/plots with the `-kapsh` name:

```
control_reset_dof  12  -kapsh
control_reset_value_constant  12  0.
```

## Parameters

| Record | Meaning |
|--------|---------|
| `materi_plasti_kappa_shear` | Adds the `kapsh` dof (size of the deviatoric plastic strain) to the node_dof records. |

## Notes

- Requires `materi_stress` and `materi_strain_plasti` for the plastic
  strain increment to exist.
- When `kapsh` is declared, the Mohr-Coulomb hardening-softening laws
  use it as the hardening variable (falling back to the total `kappa`
  dof only when `materi_plasti_kappa_shear` is absent).
- The dam_building corpus test uses it (17-dof reset including `-kapsh`).
