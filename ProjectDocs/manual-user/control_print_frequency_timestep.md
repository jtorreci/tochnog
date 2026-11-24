# control_print_frequency_timestep

## Description

`control_print_frequency_timestep` limits how often the
`control_print_*` records of the same control index run: instead of
printing at every time step, the prints happen each time after a number
of time steps has passed, and ALWAYS also at the end of the time
increment of the `control_timestep` record with the same index.

The record only affects the `control_print_*` records of the SAME
index. The exceptions `control_print`, `control_print_history` (and its
smoothing) and `control_print_data_versus_data` always run, at every
time step, independently of this record.

It should only be used in combination with `control_timestep` (with the
same index).

## Uso

```
control_timestep                10  0.04 0.41
control_print_gid               10  -yes
control_print_frequency_timestep 10  5
```

## Parámetros

| Parameter     | Meaning |
|---------------|---------|
| `10`          | Index of the control record (must match the `control_timestep` index). |
| `5`           | Number of time steps between prints. |

## Output

In the example above the `control_print_gid` (and every other gated
`control_print_*` of index 10) is executed at times 0.20, 0.40 and 0.41
- every 5 time steps (steps 5 and 10), and 0.41 because it is the end
of the time increment.

Note: the GNU time loop clamps the last partial step of an increment
into the increment end (with `0.04 0.41` the steps are 0.04..0.36 and
0.41; there is no step at 0.40). The manual example writes 0.40, which
requires a step exactly at 0.40 - use a separate final increment
(`control_timestep 10 0.04 0.40 0.01 0.01`) to reproduce it exactly.

## Example

```
control_timestep                10  0.04 0.40 0.01 0.01
control_print_dof               10  -yes
control_print_frequency_timestep 10  5
```

The file `dof.10` receives 3 blocks (at t = 0.20, 0.40, 0.41) instead
of one block per step (11).

## Validation

Test `freq_timestep` (validation-suite/test-2014): a 2D elastic model
with `control_timestep 22 0.04 0.40 0.01 0.01`, `control_print_dof 22
-separate_index`, `control_print_frd 22 -separate_sequential` and
`control_print_frequency_timestep 22 5` writes exactly 3 dof blocks and
3 .frd files whose 100CL lines carry the times 0.20, 0.40, 0.41 (the
exact times of the manual example, checked in build_safe.sh). The
control index 23 without the frequency record writes 11 blocks (one per
step). `control_print_history` on the same index still writes one line
per step (11), proving it is NOT gated.
