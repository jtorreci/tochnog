# control_print_gid_mesh_activate_gravity

Manual Professional 6.305, and its global form print_gid_mesh_activate_gravity (6.985).

## Description

Chooses which elements the GiD files hold in a calculation with `mesh_activate_gravity_time`
(elements activated from the bottom to the top).

- `-no` (default): only the elements that are active at the GiD print. The `.flavia.msh` is
  written once, at the first GiD print, so it holds only the elements already active then. The
  element results of the `.flavia.res` (`element_group`, `truss_vector`, `beam_vectors_*`) are
  written at every print for the elements active at that print.
- `-yes`: every element, at every print.

An element is active at the print time `t` when:

- no `mesh_activate_gravity_time` record selects it, or
- `t` is after its `mesh_activate_gravity_time_initial` (if set) and `t` has reached its start
  time of activation (`-method2`, the default: `t >= start`; `-method1`: `t > start`).

The start time of activation is the one of `mesh_activate_gravity_time`: interpolated between
`time_start` and `time_end` from the lowest coordinate of the element.

Nodes and nodal results are not filtered: every node is written. A model without
`mesh_activate_gravity_time` gets the same GiD files with either switch.

`control_print_gid_mesh_activate_gravity index` applies to the `control_print_gid` with the same
index and overrules the global `print_gid_mesh_activate_gravity`. Overview of the family:
[print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_mesh_activate_gravity  index  switch
print_gid_mesh_activate_gravity  switch
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `switch` | `-no` (default): only the active elements. `-yes`: every element. |

## Example

```
mesh_activate_gravity_time               10  0.  1.
control_timestep                         10  0.25  1.
control_print_gid                        10  -yes
control_print_gid_mesh_activate_gravity  10  -yes
```

## Differences with the Professional version

None measured: the rules above are those of the Professional 02-08-2026 (owner decision,
odd/tasks/gid.md L14). Two consequences are kept as the Professional writes them:

- With `-no`, a later print writes element results of elements that are not in the
  `.flavia.msh` (activated after the first print).
- If no element is active at the first print, the `.flavia.msh` is empty and the `.flavia.res`
  declares no Gauss points, although later prints write `element_group` blocks.

## Tests

`validation-suite/test-2026/` (checked by `scripts/build_safe.sh` in a scratch copy, element
numbers compared with the Professional 02-08-2026):

- `control_print_gid_mesh_activate_gravity1`: control `-no` overrules global `-yes`.
- `control_print_gid_mesh_activate_gravity2`: control `-yes` overrules global `-no`.
- `control_print_gid_mesh_activate_gravity3`: default with `-method1`.
- `control_print_gid_mesh_activate_gravity4`: default with `mesh_activate_gravity_time_initial`.
- `print_gid_mesh_activate_gravity1`: global `-yes`.
