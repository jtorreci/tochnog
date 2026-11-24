# control_print_gid_beam_vectors

Manual Professional 6.295 and control_print_gid_beam_vectors_normal 6.296.

## Description

With `-yes` the forces and moments of the `-beam` and `-truss_beam` elements are written as the
element vector results `beam_vectors_0` ... `beam_vectors_5` = f_k (t x n) on the two Gauss
points, with t the unit element axis, n the unit plotting-plane normal of
`control_print_gid_beam_vectors_normal` (default `0 0 1`) and f_k the component k of
`element_beam_force_moment`: minus the value of the first node at the first point, the value of
the second node at the second point.

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_beam_vectors         index  switch
control_print_gid_beam_vectors_normal  index  normal_x normal_y normal_z
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
control_print_gid_beam_vectors         10  -yes
```

## Differences with the Professional version

None for the vectors: equal to the Professional 02-08-2026 on a bent cantilever (element 1
beam_vectors_5 (0,-4) and (0,-2)). The scalar beam_force_moment_* results of the Professional
are not written by the GNU GiD writer.

## Tests

`validation-suite/test-2026/`: `control_print_gid_beam_vectors1`, `control_print_gid_beam_vectors_normal1` (checked by `scripts/build_safe.sh` in a scratch copy).
