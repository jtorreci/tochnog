# control_print_element_method

## Description

`control_print_element_method` selects how `control_print_element`
prints the element data:

- `-middle` (default): only the average value of the element data with
  the coordinate of the middle of the element is printed, one line per
  element.
- `-node`: the nodal values with the nodal coordinates are printed for
  each element (one line per node of the element).

It works in combination with the `control_print_element` record (same
index). If the record is not specified, `-middle` is used.

## Uso

```
control_print_element         20 -element_truss_force
control_print_element_method  20 -middle
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match the `control_print_element` index. |
| method    | `-middle` (default): average value + middle coordinate per element. `-node`: nodal values + nodal coordinates per element. |

## Output

The files written by `control_print_element` (e.g.
`element_truss_force_0.<index>`, `element_beam_force_moment_0.<index>`
for the shear force and `element_beam_force_moment_1.<index>` for the
moment). With `-middle` each element contributes one line; with `-node`
each element contributes one line per node.

## Example

```
control_print_element         20 -element_truss_force
control_print_element_method  20 -middle
control_timestep              20  1.e-1 1.e-1
```

Produces `element_truss_force_0.20` with the middle coordinate and the
truss force of every element.

## Validation

Test `elmethod` (validation-suite/test-2014): two truss elements with
known initial forces. With `-middle` the file has exactly 2 lines and
the middle coordinate equals the average of the element nodal
coordinates (2*mid == x_node2 exactly); with `-node` it has 4 lines
(one per node). The line-count and coordinate-relation A/B discriminates
the method.
