# group_beam_inertia

## Description

`group_beam_inertia` (manual Professional 6.590) gives the moments of
inertia of a `truss_beam` group:

```
group_beam_inertia index Iyy Izz J
```

- `Iyy` — moment of inertia about the local y axis;
- `Izz` — moment of inertia about the local z axis (the bending axis
  of the 2D beam in the x-y plane);
- `J` — torsional constant.

The GNU beam is the 2D x-y beam bending about the local z axis, so it
consumes `Izz` (the second value). The legacy GNU single-value form
(only `Izz`) is still accepted.

The Professional beam inputs prescribe all three rotations
(`-rotx -roty -rotz`); in a 2D model the GNU declares a single
in-plane rotation `rotz` and `-rotx`/`-roty` resolve to the same
unknown (they only appear as zero constraints at clamped ends).

## Uso

```
group_type          0  -truss_beam
group_beam_young    0  1.0
group_beam_shear    0  1.0
group_beam_inertia  0  1. 1. 1.
group_beam_memory   0  -total_linear
group_truss_area    0  1.0
group_truss_elasti_young 0  1.0
group_truss_memory  0  -total_linear
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `Iyy` | Inertia about the local y axis. |
| `Izz` | Inertia about the local z axis — used by the GNU 2D beam. |
| `J` | Torsional constant. |

`group_beam_memory` accepts `-updated`, `-updated_without_rotation`
and `-total_linear` (mapped onto the linear kinematics of the GNU
beam).

## Verification

`trubea1` and `trubea4` (cantilever beam bending tests): rc=0. A/B on
the Professional binary (25-10-2023): changing only `Izz` changes the
tip deflection by 1/I (disy = 3.333e-3 with Izz=1, 3.333e-5 with
Izz=100), confirming that the 2D bending consumes the second value.
