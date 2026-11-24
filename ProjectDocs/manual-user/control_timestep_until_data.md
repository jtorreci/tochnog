# control_timestep_until_data

## Description

`control_timestep_until_data index name_0 index_0 number_0 name_1 ...`
with `control_timestep_until_minimum index minimum_0 minimum_1 ...` and
`control_timestep_until_maximum index maximum_0 maximum_1 ...` (manual
Professional 6.390-6.392): after each timestep of control index `index`,
when a data item value is below its minimum or above its maximum, the
timesteps of that index stop; the step in which the value leaves the
range is the last one, and the calculation continues with the next
control index. Up to 10 data items; both the minimum and the maximum
record are required, one value per data item (otherwise an error).

The data item is given as for `control_print_history`: a record name,
its index and a dof name or value number, e.g. `-node_dof 4 -sigyy` or
`-post_point_dof 3 -velx`.

## Usage

```
control_timestep               10 ...
control_timestep_until_data    10 -post_point_dof 3 -velx
control_timestep_until_minimum 10 -120.
control_timestep_until_maximum 10 +120.
```

Tests: `validation-suite/test-2026/control_timestep_until1` (index 0 stops,
index 1 runs), `control_timestep_until2` (index 1 stops after one step,
index 2 runs), `control_timestep_until_error1` (maximum missing).

## Differences with the Professional version

None measured (Professional 02-08-2026: eight probes with the same step
series; the missing-maximum error is given after the first step in both).
