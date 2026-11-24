# materi_displacement_relative

## Description

`materi_displacement_relative` is an initialization-part option that adds a
**relative displacement** dof (`disrx`, `disry`, `disrz`) to the node dofs.
It tracks the displacement relative to a **reference point** in the
calculation: the current displacement minus the displacement before it was
changed by a timestep change in `control_timestep` or a displacement reset
in `control_reset_dof`.

This is handy to understand the extra displacements caused by the last
timesteps (e.g. after an excavation stage or a load increment).

Requires `materi_displacement` (which requires `materi_velocity`), plus
`materi_velocity_integrated` (the changelog: "i.c.w.
materi_velocity_integrated").

## Uso

Add it to the initialization part:

```
number_of_space_dimensions 2
derivatives
materi_velocity
materi_velocity_integrated
materi_displacement
materi_displacement_relative
materi_strain_total
materi_stress
end_initia
```

The relative displacement dofs are referenced by their labels `disrx`,
`disry`, `disrz` (e.g. in `target_item`, `control_print_history`).

## Comportamiento del punto de referencia

- **Cambio de timestep**: when a `control_timestep` record defines a
  timestep different from the previous one, the relative displacement is
  reset to 0 (the reference point moves to the current position).
- **Reset de desplazamiento**: when `control_reset_dof` resets a
  displacement dof (`disx`), the relative displacement is reset to 0.

## Related

- `materi_displacement` — the total displacement dof (`disx`/`disy`/`disz`).
- `materi_velocity_integrated` — integrated velocity (required companion).
- `control_reset_dof` — displacement reset that re-synchronizes the
  reference.
