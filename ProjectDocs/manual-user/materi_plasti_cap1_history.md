# materi_plasti_cap1_history

## Description

`materi_plasti_cap1_history` (manual Professional 4.17) is an
initialization option that adds a **scalar history variable `pc`** to
the node_dof records for the cap1 plasticity models:

> "The history variable pc for the cap1 plasticity models is added to
> the node_dof records. You need to give an initial value for it in the
> node_dof records."

Its consumer is [`group_materi_plasti_cap1`](group_materi_plasti_cap1.md)
(manual 6.691): `pc` is the position of the cap along the p axis
(`p*c = pc + c*cot(phi)`), it hardens with the cap plastic volumetric
strain and stays constant during unloading/reloading.

## Uso

In the initialization part, together with the mechanical dofs:

```
materi_velocity
materi_velocity_integrated
materi_strain_total
materi_strain_plasti
materi_plasti_cap1_history
materi_stress
```

The initial value is given via `node_dof`. The record holds **one value
per unknown slot** (with `derivatives` each physical dof occupies
`1 + ndim + 1` slots, so a 3D run needs `nuknwn = npuknwn*5` values;
the dof basename is `pc` — position 18 for the initia list above):

```
node_dof -ra -from 1 -to 8 -ra
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0
```

(`pc = 100` in the primary slot of the 19th physical dof; the other
slots are the spatial/time derivatives, zero.)

## Parámetros

| Record | Meaning |
|--------|---------|
| `materi_plasti_cap1_history` | Adds the scalar hardening variable `pc` of the cap1 models to the node_dof records. |

## Física

`pc` evolves with the cap plastic volumetric strain rate (inverted from
the manual rate form):

```
pc_dot = eps_p_cv_dot * K_ref/(lambda_star/kappa_star - 1)
         * ((pc + c*cot(phi))/p_ref)^m
```

The increment of the step is `deps_p_cv = -trace(inc_epp)` (positive in
compression, clamped to >= 0), so `pc` grows during cap compression and
stays put during unloading (elastic steps have no plastic strain) and
during dilation.

## Validation

- `mcap1.dat`: the `pc` dof starts at 100, stays at 100 during the
  elastic phase, grows 0.5/step during the 20 plastic steps to
  `pc = 109.9996` (analytic discrete fixed point 110.0, exact) and
  STAYS at that value during the 5-step elastic unload (target
  `-pc 110.4 ± 2.5`).

## Notas

- The initia is **required** by `group_materi_plasti_cap1` (the data
  check fails without it; the plasti block also exits with an error).
- The dof basename is `pc` (usable in `target_item`/`post_point_dof`
  and in plotting); like `materi_plasti_kappa`, the value is clamped to
  >= 0 after each step.
