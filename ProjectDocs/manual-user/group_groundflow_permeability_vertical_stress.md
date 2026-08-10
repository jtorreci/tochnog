# group_groundflow_permeability_vertical_stress

## Description

Stress-dependent permeability for ground water flow. The permeability
coefficient depends on the vertical effective stress `sigv`:

```
kp = a / (sigv / sig0)^b
```

where `a`, `b` and `sig0` are user parameters. The result is clamped to the
given `minimum` and `maximum` values. `sigv` is the vertical EFFECTIVE stress
(vertical total stress minus the pore pressure), with the vertical direction
taken as the gravity axis.

This record is used in combination with `group_groundflow_permeability`:
the base permeability is used where the vertical stress is extremely small
(avoiding division by zero); everywhere else the stress-dependent law above
is used.

Only intended for 2D or 3D problems.

Requires `groundflow_pressure`, `materi_stress` and
`group_groundflow_permeability`.

## Usage

```
group_groundflow_permeability <element_group>  kx [ky [kz]]
group_groundflow_permeability_vertical_stress <element_group>
                        a b sig0 minimum maximum
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `a` | Permeability scale (same units as the base permeability). |
| `b` | Stress exponent of the law. |
| `sig0` | Reference stress. |
| `minimum` | Lower clamp of the computed permeability. |
| `maximum` | Upper clamp of the computed permeability. |

## Example

```
groundflow_pressure
materi_stress
...
group_type 0  -groundflow
group_groundflow_permeability 0  0.1 0.1
group_groundflow_permeability_vertical_stress 0
                        0.1   (a)
                        1.0   (b)
                        100.  (sig0)
                        0.01  (minimum)
                        1.0   (maximum)
```

## Notes

- `sigv` is interpolated to the integration point from the nodal effective
  vertical stress. The vertical direction is the gravity axis (largest
  `force_gravity` component).
- If `|sigv|` is below a small threshold (1e-10), the base
  `group_groundflow_permeability` value is kept, avoiding division by zero.
- The effective stress uses the CURRENT pore pressure of the flow solution.
