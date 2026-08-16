# control_print_interface_stress

## Description

`control_print_interface_stress` prints the interface stresses through a
set of interfaces. In 2D the interfaces are cut by a straight line from
`(xstart,ystart)` to `(xend,yend)`, specified by
`control_print_interface_stress_2d_coordinates`. In 3D the average
interface stresses are printed at the middle of each interface element.

The stresses are written to the file `interface_stress.<index>`:
- **2D**, a line per node of each interface element:
  - first column: distance from the start point (projected on the cut
    direction)
  - second column: `interface_sign` (normal stress)
  - third column: `interface_sigt` (tangential stress)
- **3D**, a line per interface element:
  - columns 1-3: coordinates of the element middle (x y z)
  - column 4: `interface_sign` (normal stress)
  - column 5: `interface_sigt1` (first tangential stress)
  - column 6: `interface_sigt2` (second tangential stress, 3D)

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
- `control_print_interface_stress_3d_geometry name index` — 3D filter:
  only interface elements on the given geometry are printed.
- `control_print_interface_stress_3d_order order` — 3D ordering
  (`-x`, `-y`, `-z`; default element order).

## Estado de implementación

- **Implementado**: 2D — `control_print_interface_stress` +
  `control_print_interface_stress_2d_coordinates`. The normal stress
  (`interface_sign`) comes from the accumulated normal strain
  (`kn * strain_normal`); the tangential stress (`interface_sigt`) comes
  from the accumulated total tangential force
  (`ELEMENT_INTERFACE_FORCE_TANG`), so both are total accumulated
  stresses of the last converged step.
- **Implementado**: 3D — `control_print_interface_stress` + `_3d_geometry`
  + `_3d_order`. Prints the element middle, `sign`, `sigt1` and `sigt2`
  (the second tangential force comes from `ELEMENT_INTERFACE_FORCE_TANG2`).
  Validated with `iface_3d_stress` (single interface: centroid 1.5 0.5 0.5,
  sign grows with compression, sigt2 ~0) and `iface_3d_order` (two
  interfaces ordered by x).
