# group_interface_groundflow_total_pressure_tension

## Description

Interface counterpart of `group_groundflow_total_pressure_tension`. Controls
that the water pressure in an interface element is at least the value
determined from the specified `water_height`. More precisely, if the static
water pore pressure determined from the water density, the gravity and the
`water_height` exceeds the pore water pressure from the groundflow equation (in
absolute terms), this static water pressure is used instead.

This is only done if the interface normal strain (displacement difference
between the interface sides) exceeds `strain_normal_minimum`.

This option comes in handy to take care that in cracks in concrete actually the
largest water pressure from an environment is used. It so ensures that a
critical safety analysis for concrete cracking is obtained.

## Usage

```
group_interface_groundflow_total_pressure_tension <element_group>
    <strain_normal_minimum> <water_height>
```

## Parameters

| Parameter | Meaning                                                     |
|-----------|-------------------------------------------------------------|
| `element_group` | Element group, see `element_group`.                 |
| `strain_normal_minimum` | Interface normal strain above which the correction is active (crack open). |
| `water_height` | Height of the water column; the static pressure is `density * gravity * water_height`. |

## Example

```
group_type 10  -groundflow
group_interface 10  -yes
group_interface_groundflow_total_pressure_tension 10
                        1.e-3
                        5.
```

When the interface opens (normal strain > 1.e-3) under a 5 m water column, the
crack water pressure is forced to the static value.
