# control_print_gid_spring2

Manual Professional 6.313 and print_gid_spring2 6.987.

## Description

How the `-spring2` elements are written to the GiD mesh: `2` as a line with its two nodes
(`MESH Spring2Gr<group> ... ElemType Linear Nnode 2`, line `element node_1 node_2 group`), `1` as
a point on its first node (`ElemType Point Nnode 1`, line `element node_1 group`). With `2` and
[control_print_gid_group](control_print_gid_group.md) `-yes` the springs also get their
`element_group` result.

The global `print_gid_spring2` is the default of every `control_print_gid_spring2`.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_spring2  index  number_of_nodes
print_gid_spring2  number_of_nodes
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `number_of_nodes` | `1` or `2`. |

## Example

```
control_print_gid          10  -yes
control_print_gid_spring2  10  2
```

## Differences with the Professional version

- Without either record the GNU writes no spring2 elements to the GiD mesh (GiD files
  unchanged); the manual and the Professional 02-08-2026 write them with 1 node by default.

## Tests

`validation-suite/test-2026/`: `control_print_gid_spring2_1`, `control_print_gid_spring2_2`, `print_gid_spring2_1` (checked by `scripts/build_safe.sh` in a scratch copy).
