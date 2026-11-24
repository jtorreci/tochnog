# groundflow_pressure_gradient

## Description

`groundflow_pressure_gradient` (manual Professional 4.7) adds the
gradient of the hydraulic pressure head to the `node_dof` records:

```
groundflow_pressure_gradient
```

One vector dof per node, named `pres_gradx`, `pres_grady` (and
`pres_gradz` in 3D). Declare it together with
`groundflow_pressure`.

## Uso

In the initialization part:

```
number_of_space_dimensions 2
groundflow_pressure
groundflow_pressure_gradient
groundflow_saturation
groundflow_velocity
end_initia
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `pres_gradx` | dh/dx of the hydraulic pressure head. |
| `pres_grady` | dh/dy of the hydraulic pressure head. |
| `pres_gradz` | dh/dz of the hydraulic pressure head (3D). |

## Verification

`tutorial_2` parses past the initia with this record. NOTE: the
tutorial is a full dam seepage analysis whose input `include`s
`mesh.gid/mesh.dat`; the corpus harness does not stage auxiliary files
and the GNU stops at the include, so the tutorial cannot reach rc=0
under the current harness (same limitation as the earthquake tests).
