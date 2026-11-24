# control_print_gid_element_group

Manual Professional 6.301.

## Description

Only the elements of the listed element groups are written to the mesh (and their element
results). The coordinates of all nodes are still written, so the nodal results are complete.
Without the record every group is written.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_element_group  index  group_0 group_1 ...
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `group_i` | Element group number. |

## Example

```
control_print_gid                 10  -yes
control_print_gid_element_group   10  1
```

## Differences with the Professional version

None: same selection as the Professional 02-08-2026.

## Tests

`validation-suite/test-2026/`: `control_print_gid_element_group1` (checked by `scripts/build_safe.sh` in a scratch copy).
