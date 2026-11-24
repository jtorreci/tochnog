# group_materi_viscosity (and heatgeneration / user)

## Description

Newtonian (linear) viscosity model for fluids and creeping solids. The viscous
stress is

```
sigma_visc = 2*nu*D
```

where `nu` is the dynamic viscosity and `D` the rate-of-deformation tensor
(symmetric gradient of the velocity). The viscous stress is added on top of the
elastic/deviatore contribution computed by the selected material model.

Optional companions:

- `group_materi_viscosity_heatgeneration -yes` — computes the mechanical
  dissipation `2*nu*(D:D)` (heat generation rate per unit volume), used by
  coupled thermal problems.
- `group_materi_viscosity_user -yes` — delegates the viscosity to a user
  routine `user_viscosity()` (a stub in `user.cc`), so `nu` can depend on
  temperature, strain, time, etc.

Useful for ground water flow, viscous geosynthetic interfaces, creep of
asphalt/bitumen, and fluid-structure coupling in the same run.

Requires `materi_stress`, a velocity formulation (`materi_velocity`) and
`group_materi_elasti_compressibility` for nearly-incompressible flow problems.

## Usage

```
group_materi_viscosity <element_group>  nu
group_materi_viscosity_heatgeneration <element_group>  -yes|-no
group_materi_viscosity_user <element_group>  -yes|-no
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `nu` | Dynamic viscosity (stress*time units). |
| `group_materi_viscosity_heatgeneration` | `-yes` computes the heat generation `2*nu*(D:D)`. |
| `group_materi_viscosity_user` | `-yes` overrides `nu` via `user_viscosity()`. |

## Example

```
group_materi_viscosity 0  1.0
group_materi_viscosity_heatgeneration 0 -yes
```

This is the `viscos1.dat` regression test (viscosity with heat generation).

## Notes

- The user-viscosity path requires programming `user_viscosity()` in `user.cc`;
  the shipped body only prints "routine user_viscosity not programmed" and
  exits.
