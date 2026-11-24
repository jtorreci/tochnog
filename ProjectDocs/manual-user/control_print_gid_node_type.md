# control_print_gid_node_type

Manual Professional 6.307 and print_gid_node_type 6.986.

## Description

Which node coordinates are written to the GiD mesh: `-node` (the stored coordinates),
`-node_start_refined` (default, the coordinates at the start of the refinement) or
`-node_deformed_mesh` (the stored coordinates plus the displacement at the first GiD print). The
results are not changed.

The global `print_gid_node_type` is the default of every `control_print_gid_node_type`.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_node_type  index  node_type
print_gid_node_type  node_type
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `node_type` | `-node`, `-node_start_refined` (default) or `-node_deformed_mesh`. |

## Example

```
control_print_gid            10  -yes
control_print_gid_node_type  10  -node_deformed_mesh
```

## Differences with the Professional version

None: the deformed coordinates equal the Professional 02-08-2026 ones (node 2 of the test:
1.001361713241 0.004120553669).

## Tests

`validation-suite/test-2026/`: `control_print_gid_node_type1`, `print_gid_node_type1` (checked by `scripts/build_safe.sh` in a scratch copy).
