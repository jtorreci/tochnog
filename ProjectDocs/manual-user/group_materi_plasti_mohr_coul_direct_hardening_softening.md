# group_materi_plasti_mohr_coul_direct_hardening_softening

## Description

`group_materi_plasti_mohr_coul_direct_hardening_softening` (manual
Professional 6.731) is the "direct" variant of
`group_materi_plasti_mohr_coul_hardening_softening`: the same seven
parameters, but the **angles are given in DEGREES**, like every
`group_materi_plasti_mohr_coul_direct` record (the non-direct variant
uses radians).

The hardening-softening law interpolates the strength parameters
linearly with the accumulated shear plastic strain:

| Param | Meaning |
|-------|---------|
| phi_0      | friction angle at kappa_shear = 0 (degrees) |
| c_0        | cohesion at kappa_shear = 0 |
| phiflow_0  | dilatancy angle at kappa_shear = 0 (degrees) |
| phi_1      | friction angle at kappa_shear = kappa_shear_crit (degrees) |
| c_1        | cohesion at kappa_shear = kappa_shear_crit |
| phiflow_1  | dilatancy angle at kappa_shear = kappa_shear_crit (degrees) |
| kappa_shear_crit | the shear plastic strain at which the final values are reached |

## Usage

```
group_materi_plasti_mohr_coul_direct_hardening_softening 4  0. 16.5 0. 0. 12. 0. 0.18
```

with the hardening variable `kapsh` declared in the initialization part
(see `materi_plasti_kappa_shear`).

## Notes

- Same yield surface as Mohr-Coulomb; phi and c (yield AND flow) vary
  linearly with `ratio = kappa_shear / kappa_shear_crit` clamped to
  [0,1].
- The dam_building corpus test uses this record (angles 0 deg, cohesion
  16.5 -> 12, kappa_shear_crit 0.18).
