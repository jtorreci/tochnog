# materi_stress_pressure_history

## Description

`materi_stress_pressure_history` (manual Professional 4.50) is an
initialization option that adds a scalar dof to the node_dof records
containing the **maximum of the absolute value of the pressure over
time**: every step the stored value is updated to
`max(previous value, |p_new|)` with `p = -sig_mean` (positive in
compression). The value is kept while the pressure is below the maximum
(unloading/reloading) and updated when the current pressure is the new
maximum.

Its consumer is
[`group_materi_elasti_stress_pressure_history_factor`]
(group_materi_elasti_stress_pressure_history_factor.md) (manual 6.655),
which multiplies the elastic stiffness by a factor while unloading.

## Uso

In the initialization part, together with the mechanical dofs:

```
materi_velocity
materi_velocity_integrated
materi_strain_total
materi_stress
materi_stress_pressure_history
```

## Parámetros

| Record | Meaning |
|--------|---------|
| `materi_stress_pressure_history` | Adds the running maximum of `|pressure|` to the node_dof records. |

## Física

The pressure is computed from the resolved stress dofs:
`p = -(sigxx + sigyy + sigzz)/3` (positive in compression — the same
convention as the young_power/poisson_power records). The history
maximum is a true running maximum over time: it never decreases, so the
stored value marks the peak pressure ever reached at each point, which
is exactly the quantity the unloading/reloading factor needs.

## Validation

- `msph.dat` / `msph_flat.dat`: load/unload oedometer; the `sph` dof
  reaches `0.3333` (the peak `|p|`) after the loading phase and STAYS
  at `0.3333` during the unloading phase (target `-sph 0.3333` in both
  tests).

## Notas

- Requires `materi_stress` in the initialization part (the pressure is
  computed from the stress dofs).
- The dof basename is `sph` (usable in `target_item`/`post_point_dof`
  and in plotting).
