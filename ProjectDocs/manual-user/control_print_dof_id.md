# control_print_dof_id

## Description

`control_print_dof_id` adds the node number ('identity') to the files
written by `control_print_dof`. With the switch `-yes` (the default)
every line of the `dof.<index>` files contains the coordinates, the dof
value and the node number to which it belongs — e.g. in 3D:
`x y z <dof> <node_number>`. With `-no` the node number is omitted and
the classic `x y z <dof>` format is kept.

It works in combination with the `control_print_dof` record (same
index).

## Uso

```
control_print_dof    20 -separate_index
control_print_dof_id 20 -yes
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match the `control_print_dof` index. |
| switch    | `-yes` (default): also write the node number. `-no`: keep the format without the node number. |

## Output

`dof.<index>` — one line per node per dof component: `x y z <dof>
<node>` (default) or `x y z <dof>` with `-no`. The node number is the
identity of the node in the current mesh (the input node number), not a
compacted position.

## Example

```
control_print_dof               20  -separate_index
control_timestep                20  0.001 0.04
```

Produces `dof.20` with lines like `0 0 0.1152 2` (x, y, dof value,
node number) for every nodal dof component.

## Validation

Test `dof1` (validation-suite/test-2014): with the default `-yes` every
line of `dof.20` has 4 fields (x y dof node). Test `dofid_no`: with
`-no` every line has 3 fields — the A/B between the two files
discriminates the default.
