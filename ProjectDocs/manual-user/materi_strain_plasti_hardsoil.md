# materi_strain_plasti_hardsoil

## Description

`materi_strain_plasti_hardsoil` (manual Professional 4.40) adds the
plastic strain `eps_kl_plas` **specifically for the hardsoil model** to
the `node_dof` records. It is the same 6-component tensor layout as
[`materi_strain_plasti`](materi_strain_plasti.md) but with its own
dedicated dof (basenames `epphsxx epphsxy epphsxz epphsyy epphsyz
epphszz`), so it can be printed/plotted independently of the plastic
strain of other laws. The dof accumulates the plastic strain increment
of the hardsoil model (the same integration as `materi_strain_plasti`).

The manual prescribes this initialization together with
`materi_plasti_hardsoil_history` when the hardsoil model is used.

## Uso

In the initialization part:

```
materi_plasti_hardsoil_history
materi_strain_plasti_hardsoil
```

## Parámetros

None (switch-only initialization; the 6 components follow the standard
stress/strain ordering: xx xy xz yy yz zz).

## Validation

- The dof is filled by the plastic strain increments in the
  `mhardsoil_gp0_off` and `mhardsoil_plast` relaxation tests (the
  return of the hardsoil yield function) and stays zero in the purely
  elastic tests (`mhardsoil_elast`, `mhardsoil_unload`,
  `mhardsoil_gp0` with the gammap_initial control).

## Groups without the model (C107, 2026-11-04)

The dof is filled only in element groups with `group_materi_plasti_hardsoil`. In every other
group it stays 0, as in the Professional 02-08-2026. Where it is filled it
uses the nodal recovery of `materi_strain_plasti`, so with a single model
it equals `epp*` at every node. Test `materi_strain_plasti_model5`
(test-2026: druckprag on tria3, `eppdrp*` = `epp*`, `eppcap*` and
`epphs*` 0).
