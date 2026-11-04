# materi_strain_plasti_cap

## Description

`materi_strain_plasti_cap` (manual Professional 4.35) adds the plastic
strain `eps_kl_plas` **specifically for the cap models** to the
`node_dof` records. It is the same 6-component tensor layout as
[`materi_strain_plasti`](materi_strain_plasti.md) but with its own
dedicated dof (basenames `eppcaxx eppcaxy eppcaxz eppcayy
eppcayz eppcazz`), so it can be printed/plotted independently of the
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
  `eppcazz = -0.036` EXACT (analytic: deps_p_cv = 0.0054/step x 20 /
  3 = 0.036; measured -0.0360001). The elastic twin
  `mstrain_cap_elast` gives exactly 0. A/B: without the initia the
  dof does not exist (target `-eppcazz` fails).

## Groups without the model (C107, 2026-11-04)

The dof is filled only in element groups with `group_materi_plasti_cap` or `group_materi_plasti_cap1`. In every other
group it stays 0, as in the Professional 02-08-2026. Where it is filled it
uses the nodal recovery of `materi_strain_plasti`, so with a single model
it equals `epp*` at every node. Test `materi_strain_plasti_model5`
(test-2026: druckprag on tria3, `eppdr*` = `epp*`, `eppca*` and
`eppha*` 0).
