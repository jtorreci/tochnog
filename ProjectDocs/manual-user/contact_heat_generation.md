# contact_heat_generation

## Description

Factor determining how much of the frictional energy loss at a slipping
contact is transformed into heat on the temperature degrees of freedom:
`Q = eta * Ff * vf`, where `Ff` is the friction force, `vf` the slip
velocity and `eta` this factor (manual Professional, sections 2.3.2 and
6.100). The heat is injected into the `condif_temperature` equation of
the contacter node (full amount when the target is an analytic
`contact_geometry` face, half to contacter and half to the target
otherwise).

This only makes sense when friction is active (`contact_friction` or
`contact_plasti_friction`) and `condif_temperature` is initialized.
The factor should be between 0 and 1. Default 0 (no heat).

`contact_heat_generation` is the Professional keyword name. The legacy
GNU name `contact_heatgeneration` (without underscore) remains accepted;
if both records exist, the Professional name wins.

## Usage

```
contact_heat_generation <factor>
```

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `factor`  | Fraction of frictional energy converted to heat. Default `0`.  |

## Example

```
contact_geometry 0  -geometry_line 1
contact_penalty_velocity 100000.
contact_friction 0.5
contact_heat_generation 0.5
```

A block sliding on the contact face while pressed against it converts
half of the frictional loss into heat. Verified via debug trace: the
internal `friction_energy` scales exactly with the factor (0.5 -> 1.312,
1.0 -> 2.624 in the reference probe).
