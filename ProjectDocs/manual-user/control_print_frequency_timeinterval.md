# control_print_frequency_timeinterval

## Description

`control_print_frequency_timeinterval` limits how often the
`control_print_*` records of the same control index run: instead of
printing at every time step, the prints happen each time after a time
interval has passed, and ALWAYS also at the end of the time increment
of the `control_timestep` record with the same index.

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
control_print_frequency_timeinterval 10  0.15
```

## Parámetros

| Parameter     | Meaning |
|---------------|---------|
| `10`          | Index of the control record (must match the `control_timestep` index). |
| `0.15`        | Time interval in the same time units as `control_timestep`. |

## Output

In the example above the `control_print_gid` (and every other gated
`control_print_*` of index 10) is executed at times 0.16, 0.32 and 0.41
- the first time 0.15 time units have passed since the start of the
increment, the second after another 0.15, and 0.41 because it is the
end of the time increment (10 steps of 0.04 plus the last partial step
clamped to the increment end).

## Example

```
control_timestep                10  0.04 0.41
control_print_dof               10  -yes
control_print_frequency_timeinterval 10  0.15
```

The file `dof.10` receives 3 blocks (at t = 0.16, 0.32, 0.41) instead
of one block per step (10).

## Validation

Test `freq_timeint` (validation-suite/test-2014): a 2D elastic model
with `control_timestep 10 0.04 0.41`, `control_print_dof 10
-separate_sequential`, `control_print_frd 10 -separate_sequential` and
`control_print_frequency_timeinterval 10 0.15` writes exactly 3 dof
files and 3 .frd files whose 100CL lines carry the times 0.16, 0.32,
0.41 (the exact times of the manual example, checked in build_safe.sh).
The control index 11 without the frequency record writes 10 files (one
per step). `control_print_history` on the same index still writes one
line per step (10), proving it is NOT gated.
