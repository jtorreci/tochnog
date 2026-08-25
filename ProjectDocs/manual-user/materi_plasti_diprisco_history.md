# materi_plasti_diprisco_history

## Description

`materi_plasti_diprisco_history number_of_history_variables` (manual
Professional 4.18) adds the history variables of the di Prisco
plasticity models to the `node_dof` records. It is the per-model name
of [`materi_history_variables`](materi_history_variables.md): the same
mechanism, the same shared history dof (basenames `hisv0 .. hisvN-1`).
The manual prescribes:

- `materi_plasti_diprisco_history 11` for
  `group_materi_plasti_diprisco` (the chi tensor 9 values + beta + rc);
- `materi_plasti_diprisco_history 12` for
  `group_materi_plasti_diprisco_density` (adds the relative density).

## Uso

In the initialization part:

```
materi_plasti_diprisco_history 11
```

The initial history values are given via `node_dof` (e.g. for a
normally consolidated sand: chi = delta_ij/sqrt(3), beta = 0.0001 and
rc = sqrt(3) times the mean pressure).

## Parámetros

| # | Meaning |
|---|---------|
| number_of_history_variables | number of history variables (11 for diprisco, 12 for diprisco_density) |

## Validation

- `mdiprisco_hist` (diprisc1 rig with the alias initia instead of
  `materi_history_variables 11`): the target `hisv10 = -152.795`
  matches diprisc1 EXACTLY (same physics, same history dof).
- `mstrain_diprisco` also uses the alias and keeps the diprisco
  history target `hisv10 = -152.795`.
