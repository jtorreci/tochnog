# control_print_number_iterations

## Description

`control_print_number_iterations` is a console monitor: with the switch
`-yes` Tochnog prints the iteration number on stdout during the
equilibrium iterations of every time step. This is convenient in very
large calculations, where you want to monitor the evolution of the
calculation without opening result files.

It is NOT the data item `-inverse_iteration_number` (which stores the
number of iterations of the already finished step for the inverse
analysis): the monitor prints the live counter during the iterations.

## Uso

```
control_print_number_iterations 20 -yes
control_timestep_iterations     20 8
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. |
| switch    | `-yes`: print the iteration counter during the equilibrium iterations. |

## Output

One line per equilibrium iteration on stdout:

```
control_print_number_iterations: time <t> iteration <n>
```

The lines are a monitor only — no result file is written and no record
is stored.

## Example

```
control_print_number_iterations 20 -yes
control_timestep_iterations     20 8
control_timestep                20  1.e-1 2.e-1
```

With 2 time steps of 8 iterations each, 16 monitor lines are printed
(the `numit` test).

## Validation

Test `numit` (validation-suite/test-2014): the captured stdout of the
run contains exactly 16 `control_print_number_iterations:` lines
(2 steps x 8 iterations — the count discriminates both the switch and
the per-iteration printing).
