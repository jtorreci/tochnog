# control_print_interface_stress

## Description

`control_print_interface_stress` prints the interface stresses through a
set of interfaces. In 2D the interfaces are cut by a straight line from
`(xstart,ystart)` to `(xend,yend)`, specified by
`control_print_interface_stress_2d_coordinates`.

The stresses are written to the file `interface_stress.<index>`:
- first column: distance from the start point (projected on the cut
  direction)
- second column: `interface_sign` (normal stress)
- third column: `interface_sigt` (tangential stress)

A line is written for each node of each interface element.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_timestep` record:

```
control_print_interface_stress               0  -separate_index
control_print_interface_stress_2d_coordinates 0  0. 0. 3. 0.
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `0`       | Index of the control record. Must match the `control_timestep` index. |
| switch    | `-separate_index`: file `interface_stress.<index>`. `-separate_sequential`: files `interface_stress.0`, `interface_stress.1`, ... |

## Related

- `control_print_interface_stress_2d_coordinates index xstart ystart
  xend yend` — cut line (2D only).
- `control_print_interface_stress_3d_geometry`, `_3d_order` — 3D variants
  (not implemented yet).

## Estado de implementación

- **Implementado**: 2D — `control_print_interface_stress` +
  `control_print_interface_stress_2d_coordinates`. The normal stress
  (`interface_sign`) comes from the accumulated normal strain
  (`kn * strain_normal`); the tangential stress (`interface_sigt`) comes
  from the accumulated total tangential force
  (`ELEMENT_INTERFACE_FORCE_TANG`), so both are total accumulated
  stresses of the last converged step.
- **Pendiente**: 3D (`_3d_geometry`, `_3d_order`).
