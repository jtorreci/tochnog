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
