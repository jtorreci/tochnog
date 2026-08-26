# control_print_dof_rhside

## Description

`control_print_dof_rhside` prints the right-hand-side of the primary
dofs together with the coordinates at which they hold. For example the
file `temp_rhside.<index>` contains lines with `x y z` and the
right-hand-side of `-temp` (that is, the heat flux); in 1D only `x` is
printed, etc.

This is the Professional name of the GNU keyword
`control_print_unknownsrhside` — both names activate the same output.

## Uso

```
control_print_dof_rhside 20 -yes
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. |
| switch    | `-yes`: print the right-hand-side files. `-no`/omitted: no output. |

## Output

One file per primary dof, named `<dof>_rhside.<index>` (e.g.
`temp_rhside.20`, `velx_rhside.20`), with one `x y z <rhs>` line per
node. The right-hand-side is the residual of the corresponding equation.

## Example

```
control_print_dof_rhside 20 -yes
control_timestep        20  0.1 0.1
```

Produces `velx_rhside.20` / `vely_rhside.20` with the force balance
per node.

## Validation

Test `dofrhside` (validation-suite/test-2014): smoke — the file
`velx_rhside.20` is generated and parseable (x y rhs columns in 2D).
