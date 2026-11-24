# control_print_dof_smooth_dof

## Description

`control_print_dof_smooth_dof` smooths the results printed by
[`control_print_dof`](control_print_dof.md) with the same index (manual
Professional 6.271). With `dof_0 dof_1 ...` you specify the dofs to be
smoothed (dof labels like `-velx`, `-disy`); the special option `-all`
smooths every dof.

The smoothing is NODAL: each pass replaces the value of a node by the
average of the values of its NEIGHBOUR nodes (nodes connected to it by
at least one element). The smoothing is applied to the values written to
the `dof.<index>` file, BEFORE the lines are written; the node id column
([`control_print_dof_id`](control_print_dof_id.md)) is not smoothed.

The smoothing is done a number of times, with increasingly smooth
results. The number of passes comes from
[`control_print_dof_smooth_n`](control_print_dof_smooth_n.md); if it is
not given, 10 passes are done.

## Uso

Place it in the data part with the same `icontrol` as the
`control_print_dof` record:

```
control_print_dof            0  -separate_index
control_print_dof_smooth_dof 0  -all
control_print_dof_smooth_n   0  3
control_timestep            0  1.0 1.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `0`       | Index of the control record. Must match `control_print_dof`. |
| `dof_0 ...` | Dof labels to smooth (e.g. `-velx -vely`). `-all`: every dof. |

## Output

No output by itself; it modifies the values written to `dof.<index>` by
[`control_print_dof`](control_print_dof.md).

## Example

1D bar chain with the velocities prescribed to `0, 1, 2, 3, 4` at
`x = 0..4`; one smoothing pass gives each interior node the exact
average of its two neighbours (`1, 1, 2, 3, 3`), the end nodes take
their single neighbour:

```
control_print_dof            60  -separate_index
control_print_dof_id         60  -no
control_print_dof_smooth_dof 60  -all
control_print_dof_smooth_n   60  1
control_timestep            60  1.0 1.0
```

Produces `dof.60` with `x <value>` lines `0 1 / 1 1 / 2 2 / 3 3 / 4 3`.

## Differences with the Professional version

- The node itself is NOT included in the average (only the neighbour
  nodes). A node without elements (no neighbours) keeps its value.
- Dof labels that do not exist in the model are skipped.
