# group_materi_expansion_volume

## What it does

`group_materi_expansion_volume index beta` (manual Professional 6.665) is
the volume expansion coefficient of a `-materi` element group in a model
with `condif_temperature`. The material density becomes

    rho = rho0 * ( 1 - beta * T )

with `rho0` the `group_materi_density` and `T` the temperature dof. The
reduced density is used for the gravity force AND, since 2026-10-06
(audit K04, Pro test `incnav3`), for the inertia of the momentum equation
(the mass), as the Professional 02-08-2026.

## Syntax

```
group_materi_expansion_volume  index  beta
```

- `index`: element group.
- `beta`: volume expansion coefficient (1/temperature).

## Example (Pro test incnav3)

```
condif_temperature
materi_velocity
...
group_materi_density 0  1.0
group_materi_expansion_volume 0  1.
force_gravity 0. 1.
```

With T = 0.5 the density is 0.5; the expansion and gravity forces give a
velocity rate of -2.5/0.5 = -5 (`tvely` -4.9904 after one step of 0.01).
Before 2026-10-06 the GNU used the mass of rho0 (-2.4976). Test:
`validation-suite/test-2026/expansion_volume_inertia1.dat`.

## Differences with the Professional version

None known. The legacy sfnet `incnav3` target (-2.5) predates this
convention (BACKLOG B11).
