# control_print_gid_coord

Manual Professional 6.298 and print_gid_coord 6.982.

## Description

With `-yes` the GiD results get the nodal vector result `coord` (the node coordinates) at
every print. With `-no` it is not written.

The global `print_gid_coord` is the default of every `control_print_gid_coord`; a
`control_print_gid_coord` overrules it for the `control_print_gid` with the same index.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_coord  index  switch
print_gid_coord  switch
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `switch` | `-yes`: write the `coord` result; `-no`: do not. |

## Example

```
control_print_gid        10  -yes
control_print_gid_coord  10  -yes
```

## Differences with the Professional version

- Without either record the GNU writes no `coord` result: the GiD files stay the ones the GNU
  always wrote. The manual and the Professional 02-08-2026 write it by default (`-yes`).
- Format of the GNU GiD writer (`Result "coord" "Load Analysis" <time> Vector OnNodes`).

## Tests

`validation-suite/test-2026/`: `control_print_gid_coord1`, `control_print_gid_coord2`, `print_gid_coord1` (checked by `scripts/build_safe.sh` in a scratch copy).
