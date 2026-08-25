# materi_strain_plasti_druckprag

## Description

`materi_strain_plasti_druckprag` (manual Professional 4.39) adds the
plastic strain `eps_kl_plas` **specifically for the drucker-prager
model** (`group_materi_plasti_druck_prag`) to the `node_dof` records.
It is the same 6-component tensor layout as
[`materi_strain_plasti`](materi_strain_plasti.md) but with its own
dedicated dof (basenames `eppdrpxx eppdrpxy eppdrpxz eppdrpyy eppdrpyz
eppdrpzz`). The dof accumulates the plastic strain increment of the
drucker-prager law (the same integration as `materi_strain_plasti`).

## Uso

In the initialization part:

```
materi_strain_plasti_druckprag
```

## Parámetros

None (switch-only initialization; the 6 components follow the standard
stress/strain ordering: xx xy xz yy yz zz).

## Validation

- `mstrain_druckprag` (pure shear rig of mdp_shear, phi = 0):
  `eppdrpxy = 0.910768` IDENTICAL to the generic `materi_strain_plasti`
  eppxy of the same run (the dedicated dof records the same plastic
  shear strain as the generic tensor); the stress target
  `sigxy = 34.64` (2c/sqrt(3)) still holds. The elastic twin
  `mstrain_druckprag_elast` (c = 1e6) gives exactly 0. A/B: without
  the initia the dof does not exist.
