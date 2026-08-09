# bounda_water

## Description

Applies the hydrostatic pore-water pressure as a prescribed boundary condition.
When a dof is prescribed with value `-pres` (pore pressure) and `bounda_water`
is active, the pressure is not taken from a fixed value but computed from the
height of the water column between the node and the phreatic level:

```
pres = density_water * g * (water_level - y)
```

where `g` is the vertical component of the gravity vector and `y` is the
vertical coordinate of the node. The result is a linear hydrostatic pressure
profile below the phreatic level.

Requires the `groundflow_pressure` and `groundflow_velocity` models and the
`bounda_dof`/`bounda_unknown` record prescribing the `-pres` dof on the same
nodes.

## Usage

```
bounda_water <index> <yes|no>
```

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `index`   | Record index; matches the index of the `bounda_unknown`/`bounda_dof` record. |
| `yes|no`  | `-yes` activates the hydrostatic water pressure; `-no` disables it. |

The pressure depends on the following (non-indexed) data records:

| Record                        | Meaning                                             |
|-------------------------------|-----------------------------------------------------|
| `groundflow_density`          | Density of the water.                               |
| `force_gravity`               | Gravity vector; only the vertical component is used.|
| `groundflow_phreaticlevel`    | Phreatic (water) level; the first value is used.    |
| `bounda_time`                 | Multiplier of the computed static pressure.         |

## Example

```
groundflow_density 1000.
force_gravity 0. -9.81 0.
groundflow_phreaticlevel 5.
bounda_time 1. -1.
bounda_unknown 1 pres pres pres pres
bounda_water 0 -yes
```

Nodes below the phreatic level (`y < 5.`) get the pore pressure
`1000. * 9.81 * (5. - y)`, scaled by `bounda_time`.

## Note on the pressure value

`bounda_water` applies the **full computed static pressure** (including
positive/compression values below the phreatic level). It does NOT cap the
pressure at `groundflow_pressure_atmospheric`.

This is intentional, but note that the internal function
`groundflow_phreatic_coord()` (used elsewhere by the groundflow solver) caps the
static pressure at `groundflow_pressure_atmospheric`, which defaults to **0**.
With that default, that function only returns meaningful pressures for
**suction** (negative pressures, nodes above the phreatic level). If your model
expects the capped behaviour, define `groundflow_pressure_atmospheric`
explicitly so both paths agree.
