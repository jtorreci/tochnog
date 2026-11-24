# control_print_gid_truss_vector

Manual Professional 6.314 and control_print_gid_truss_vector_normal 6.315.

## Description

With `-yes` the normal force N of the `-truss` and `-truss_beam` elements is written as the
element vector result `truss_vector` = N (t x n) on the two Gauss points of the element, with t
the unit element axis (first to second node) and n the unit plotting-plane normal of
`control_print_gid_truss_vector_normal` (default `0 0 1`). A normal parallel to the element gives
zero vectors. The scalar truss force is not written by the GNU GiD writer.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_truss_vector         index  switch
control_print_gid_truss_vector_normal  index  normal_x normal_y normal_z
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `switch` | `-yes` or `-no` (default). |
| `normal_*` | Plotting-plane normal, default `0 0 1`. |

## Example

```
control_print_gid                      10  -yes
control_print_gid_truss_vector         10  -yes
control_print_gid_truss_vector_normal  10  0 0 -1
```

## Differences with the Professional version

None for the vector: the values equal the Professional 02-08-2026 (element at 45 degrees with
N = 1.414213562373: (1,-1)).

## Tests

`validation-suite/test-2026/`: `control_print_gid_truss_vector1`, `control_print_gid_truss_vector_normal1` (checked by `scripts/build_safe.sh` in a scratch copy).
