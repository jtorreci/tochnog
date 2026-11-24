# group_groundflow_total_pressure_tension

## Description

Controls that the water pressure in an element is at least the value
determined from the specified `water_height`. More precisely, if the static
water pore pressure determined from the water density, the gravity and the
`water_height` exceeds the pore water pressure from the groundflow equation (in
absolute terms), this static water pressure is used instead.

This is only done if the largest eigenvalue of `materi_strain_plastic_tension`
exceeds `plastic_tension_minimum`. To calculate the eigenvalues of
`materi_strain_plastic_tension` include
`post_calcul -materi_strain_plasti_tension -prival` in the input file.

This option comes in handy to take care that in cracks in concrete actually the
largest water pressure from an environment is used. It so ensures that a
critical safety analysis for concrete cracking is obtained.

## Usage

```
group_groundflow_total_pressure_tension <element_group>
    <plastic_tension_minimum> <water_height>
```

## Parameters

| Parameter | Meaning                                                     |
|-----------|-------------------------------------------------------------|
| `element_group` | Element group, see `element_group`.                 |
| `plastic_tension_minimum` | Lower bound of the largest eigenvalue of `materi_strain_plastic_tension` above which the correction is active. |
| `water_height` | Height of the water column; the static pressure is `density * gravity * (water_height - coord)`. |

## Example

```
group_groundflow_total_pressure_tension 0
                        0.001
                        5.
```

With a plastic crack (eigenvalue > 0.001) and a 5 m water column above the
domain, the crack water pressure is forced to the static value.
