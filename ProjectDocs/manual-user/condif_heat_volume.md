# condif_heat_volume

## Description

Distributed volume heat source S for the heat transfer (condif)
equation (manual Professional 6.83): contributes `S` per unit volume to
the temperature equation (rhs += ∫ S N_i dV).

Companions (same index):

- `condif_heat_volume_element` — restrict to the listed elements.
- `condif_heat_volume_element_group` — restrict to elements of the
  listed element groups.
- `condif_heat_volume_geometry` — restrict to elements with ALL nodes in
  the geometry.
- `condif_heat_volume_factor` — spatial polynomial
  `a0 + a1*x + a2*x^2 + ...` at the integration point.
- `condif_heat_volume_time` — `time load` pairs (linear interpolation).
- `condif_heat_volume_sine` — `start_time freq_0 amp_0 freq_1 amp_1...`.
- `condif_heat_volume_user` + `condif_heat_volume_user_parameters` —
  heat defined by the user routine `user_condif_heat_volume`
  (edit user.cc; enabling `-yes` without programming it is an error).

`condif_temperature` must be an unknown.

## Usage

```
condif_heat_volume <index> <heat>
```

## Example

```
condif_heat_volume 0  1.0
condif_heat_volume_element 0  2
```

1D bar (nodes 0,1,2; T=0 at both ends; k=1) heated only in element 2:
steady state gives T(center) = 0.25 (analytic); heating both elements
would give 0.5 (test `condif_heat_vol`). With
`condif_heat_volume_factor 0 0. 1.` the source is S(x)=x and
T(center) = 0.5 (test `condif_heat_vol2`).
