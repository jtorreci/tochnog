# materi_strain_plasti_cap

## Description

`materi_strain_plasti_cap` (manual Professional 4.35) adds the plastic
strain `eps_kl_plas` **specifically for the cap models** to the
`node_dof` records. It is the same 6-component tensor layout as
[`materi_strain_plasti`](materi_strain_plasti.md) but with its own
dedicated dof (basenames `eppcapxx eppcapxy eppcapxz eppcappyy
eppcapyz eppcapzz`), so it can be printed/plotted independently of the
plastic strain of other laws. The dof accumulates the plastic strain
increment of the active cap law (the same integration as
`materi_strain_plasti`).

## Uso

In the initialization part:

```
materi_plasti_cap1_history
materi_strain_plasti_cap
```

## Parámetros

None (switch-only initialization; the 6 components follow the standard
stress/strain ordering: xx xy xz yy yz zz).

## Validation

- `mstrain_cap` (isotropic compression with cap1, 20 plastic steps):
  the dedicated dof accumulates the cap plastic volume strain:
  `eppcapzz = -0.036` EXACT (analytic: deps_p_cv = 0.0054/step x 20 /
  3 = 0.036; measured -0.0360001). The elastic twin
  `mstrain_cap_elast` gives exactly 0. A/B: without the initia the
  dof does not exist (target `-eppcapzz` fails).
