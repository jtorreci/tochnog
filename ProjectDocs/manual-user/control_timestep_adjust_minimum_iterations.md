# control_timestep_adjust_minimum_iterations

## Description

`control_timestep_adjust_minimum_iterations index switch` (manual
Professional 6.384): with `-yes` (Professional default) the program may
increase the minimum number of iterations of the automatic time stepping
when it thinks that helps; `-no` keeps the minimum (2 or the
[control_timestep_iterations_automatic_minimum_maximum_wished](control_timestep_iterations_automatic_minimum_maximum_wished.md)
value).

## Usage

```
control_timestep_adjust_minimum_iterations index switch
```

## Differences with the Professional version

The record is accepted, but the GNU always keeps the minimum (as with
`-no`): the manual does not specify the adjustment rule (BACKLOG C94).
