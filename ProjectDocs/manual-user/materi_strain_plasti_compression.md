# materi_strain_plasti_compression

## Description

`materi_strain_plasti_compression` (manual Professional 4.36) adds the
plastic strain `eps_kl_plas` **specifically for the compression model**
(`group_materi_plasti_compression`) to the `node_dof` records. It is
the same 6-component tensor layout as
[`materi_strain_plasti`](materi_strain_plasti.md) but with its own
dedicated dof (basenames `eppcmpxx eppcmpxy eppcmpxz eppcmpyy eppcmpyz
eppcmpzz`). The dof accumulates the plastic strain increment of the
compression law (the same integration as `materi_strain_plasti`).

## Uso

In the initialization part:

```
materi_strain_plasti_compression
```

## Parámetros

None (switch-only initialization; the 6 components follow the standard
stress/strain ordering: xx xy xz yy yz zz).

## Validation

- `mstrain_compression` (confined compression, sig_yield = 1.0, 10
  steps): the dedicated dof fills `eppcmpzz = -8.96e-05`. The elastic
  twin `mstrain_compression_elast` (load far below yield) gives
  exactly 0. A/B: without the initia the dof does not exist.

## Groups without the model (C107, 2026-11-04)

The dof is filled only in element groups with `group_materi_plasti_compression` or `group_materi_plasti_compression_direct`. In every other
group it stays 0, as in the Professional 02-08-2026. Where it is filled it
uses the nodal recovery of `materi_strain_plasti`, so with a single model
it equals `epp*` at every node. Test `materi_strain_plasti_model5`
(test-2026: druckprag on tria3, `eppdrp*` = `epp*`, `eppcap*` and
`epphs*` 0).
