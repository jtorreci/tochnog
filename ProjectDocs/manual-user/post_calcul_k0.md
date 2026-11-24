# post_calcul -k0 operator

## Description

`post_calcul` with the `-k0` operator on the `-materi_stress` matrix
computes the ratio of the average horizontal stress over the vertical
one (manual Professional 6.901, "Specially for geotechnics"): the
earth-pressure coefficient at rest of the resolved stress state:

- 2D: `0.5*(sigxx + sigzz) / sigyy`
- 3D: `0.5*(sigxx + sigyy) / sigzz`

One scalar value per node/point is produced. The output item is named
`k0_sig` (the `post_calcul_label` of the Professional `.dbs`), so the
results are stored in the `node_dof_calcul` / `post_point_dof_calcul`
records under the slot you can select in a `target_item` with
`-k0_sig`.

## Usage

```
post_calcul -materi_stress -k0
target_item <n> -post_point_dof_calcul <point> -k0_sig
```

## Example

See `k0.dat`: a laterally confined column under gravity with
nu = 0.1 resolves `sigxx = sigzz = nu/(1-nu) * sigyy`, so
`k0_sig = 0.11111111` at every node (target, verified against the
Professional `.dbs`).
