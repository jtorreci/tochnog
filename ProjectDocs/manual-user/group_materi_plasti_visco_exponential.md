# group_materi_plasti_visco_exponential / _power / _always

## Description

Viscoplastic (rate-dependent) plasticity for soils and other geomaterials.
These records must be used in combination with a plasticity model
(`group_materi_plasti_*`, e.g. Mohr-Coulomb, von Mises, Cam-Clay). They replace
the purely incremental plastic strain with a rate law driven by the yield
function value `f` and the mean (compressive) pressure.

Two rate laws are implemented:

- **Exponential** (`group_materi_plasti_visco_exponential gamma alpha`):
  ```
  lambda = gamma * (-pressure) * exp(alpha*f)
  ```
  with the exponent argument clamped to `EPS_VISCO = 3.` (the professional
  default of `group_materi_plasti_visco_exponential_limit`).

- **Power** (`group_materi_plasti_visco_power eta p f_ref`):
  ```
  lambda = eta * (f/f_ref)^p
  ```
  (`f_ref != 0` is asserted).

`lambda` is the viscoplastic multiplier used by the flow rule. Both laws give
zero strain rate on the yield surface (`f=0`) and rate-dependent flow when
`f>0` (overstress), which regularizes strain localization and models creep of
soft soils.

Optional `group_materi_plasti_visco_always -yes` forces plastic
(viscoplastic) flow even when `f < EPS_F` (i.e. plasticity is applied at every
iteration, not only when the stress exceeds the yield surface).

Requires `materi_stress`, `materi_strain_plasti`, a plasticity model, and a
velocity formulation.

## Usage

```
group_materi_plasti_visco_exponential <element_group>  gamma alpha
group_materi_plasti_visco_power      <element_group>  eta p f_ref
group_materi_plasti_visco_always     <element_group>  -yes|-no
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `gamma` | Fluidity parameter of the exponential law (>= 0). |
| `alpha` | Exponential exponent coefficient (>= 0). |
| `eta`  | Viscoplastic viscosity coefficient of the power law (>= 0). |
| `p`    | Power exponent of the overstress ratio. |
| `f_ref` | Reference stress for the overstress ratio (must be non-zero). |
| `visco_always` | `-yes` applies viscoplastic flow at every iteration. |

`gamma<=0`, `alpha<=0` and `eta<=0` raise `db_error()` at the constitutive
level.

## Example

```
group_materi_plasti_visco_exponential 0  0.01  0.5
```

## Notes

- The exponent argument limit is hardcoded to 3 in `stress.cc`
  (`EPS_VISCO`); the professional keyword
  `group_materi_plasti_visco_exponential_limit` is NOT implemented.
- The `_name`/`_values` variants (per-plasticity-model parameters) of the
  professional manual are NOT implemented.
