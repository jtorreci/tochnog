# control_print_gid_element_mpc

Manual Professional 6.302.

## Description

With `-no` the elements that have an mpc'd node are not written to the GiD mesh (nor their
element results). An mpc'd node is the first (slave) node of an `mpc_node_number` record; an
element with only master nodes is written. Default `-yes`: every element is written.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_element_mpc  index  switch
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `switch` | `-yes` (default) or `-no`. |

## Example

```
mpc_node_number 0  7 -vely  6 -vely
control_print_gid              10  -yes
control_print_gid_element_mpc  10  -no
```

## Differences with the Professional version

None: measured on the Professional 02-08-2026 (slave node 7 drops elements 2 and 3, element 1
with the master node 6 stays).

## Tests

`validation-suite/test-2026/`: `control_print_gid_element_mpc1` (checked by `scripts/build_safe.sh` in a scratch copy).
