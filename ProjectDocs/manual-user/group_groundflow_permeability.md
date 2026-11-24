# group_groundflow_permeability

## Description

Permeability coefficient of the ground water flow in each space
direction (manual Professional 6.618):

```
group_groundflow_permeability <index> <pex> [<pey> [<pez>]]
```

- In 1D specify only `pex`, in 2D `pex pey`, in 3D `pex pey pez`.
- If you specify only ONE value, that value is used in EVERY space
  direction (isotropic permeability).

`<index>` is the element group of the groundflow material (see
`group_type ... -groundflow`).

## Usage

```
group_type 0 -materi -groundflow
group_groundflow_permeability 0 1.e-2 1.e-2      (2D, anisotropic input)
group_groundflow_permeability 0 1.01937e-4        (isotropic: k/gamma_w)
```

The record is a companion of the storage equation
`c * h_dot = kp_i * d2h/dxi2 + div(v)` (manual 2.4.1): the `kp_i` are
the `group_groundflow_permeability` values per direction. For the
stress-dependent variant see
[group_groundflow_permeability_vertical_stress](group_groundflow_permeability_vertical_stress.md).

## Notes

- Ground flow velocities, when initialized with `groundflow_velocity`,
  follow from `v_g,i = kp_i * (-d h / d x_i)` (with the sign convention
  of the manual); the isotropic value fills all `kp_i`.
- The vertical-stress law falls back to this record where the vertical
  effective stress is extremely small (avoiding division by zero).
- Units: the corpus tests give the permeability as `k/gamma_w` (hydraulic
  conductivity over water unit weight, Darcy flux per pressure gradient).
