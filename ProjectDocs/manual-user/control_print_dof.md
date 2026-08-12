# control_print_dof

## Description

`control_print_dof` prints the values of the primary dofs together with
the coordinates at which they hold. It is the low-level companion of the
format-based outputs (`control_print_tabular`, `control_print_vtk`, ...):
plain ASCII files with one `x y z <dof>` line per node, useful for
external post-processing or quick inspection.

The coordinates themselves are also written to a separate file (`coord.*`),
once (on the first call).

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_timestep` record:

```
control_print_dof 20 -separate_index
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match an active `control_timestep`. |
| switch    | `-yes`: file `dof.<index>`. `-separate_index`: file `dof.<index>`. `-separate_sequential`: files `dof.0`, `dof.1`, ... one per call. |

## Output

- `dof.<index>` — one line per node per dof component:
  `x y z <dof_value>` (only `x` in 1D, `x y` in 2D, `x y z` in 3D). Each
  call to `control_print_dof` appends the current state, so the file
  grows as a time series (one block of all dofs per time step).
- `coord.<index>` — the nodal coordinates, written only on the first call.

Vectors are printed per component (`velx`, `vely`, ...), matrices (stress)
per Voigt component (`sigxx`, `sigyy`, `sigzz`, `sigxy`, `sigxz`,
`sigyz`).

## Example

```
control_print_dof               20  -separate_index
control_timestep                20  0.001 0.04
```

Produces `dof.20` with, for each time step, one block of
`x y <dof_value>` lines for every nodal dof component, and `coord.20`
with the coordinates.
