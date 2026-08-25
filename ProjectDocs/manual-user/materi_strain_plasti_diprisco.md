# materi_strain_plasti_diprisco

## Description

`materi_strain_plasti_diprisco` (manual Professional 4.37) adds the
plastic strain `eps_kl_plas` **specifically for the di Prisco model**
(`group_materi_plasti_diprisco`) to the `node_dof` records. It is the
same 6-component tensor layout as
[`materi_strain_plasti`](materi_strain_plasti.md) but with its own
dedicated dof (basenames `eppdipxx eppdipxy eppdipxz eppdipyy eppdipyz
eppdipzz`). The dof accumulates the plastic strain increment of the
di Prisco law (the same integration as `materi_strain_plasti`).

## Uso

In the initialization part (with the di Prisco history variables, see
[`materi_plasti_diprisco_history`](materi_plasti_diprisco_history.md)):

```
materi_plasti_diprisco_history 11
materi_strain_plasti_diprisco
```

## Parámetros

None (switch-only initialization; the 6 components follow the standard
stress/strain ordering: xx xy xz yy yz zz).

## Validation

- `mstrain_diprisco` (undrained axisymmetric triaxial, diprisc1 rig):
  the dedicated dof fills `eppdipyy = -0.0105714` (axial compression)
  and `eppdipxx = eppdipzz = +0.00350735` (radial = hoop); the
  di Prisco history target `hisv10 = -152.795` matches diprisc1. The
  elastic twin `mstrain_diprisco_elast` gives exactly 0. A/B: without
  the initia the dof does not exist.
