# control_print_beam_force_moment_coordinates

## Description

Defines the cut segment for
[`control_print_beam_force_moment`](control_print_beam_force_moment.md)
(manual Professional 6.263). The forces and moments of the beams
crossed by the segment from `(xstart, ystart, zstart)` to
`(xend, yend, zend)` are printed.

## Input syntax

```
control_print_beam_force_moment_coordinates <index> xstart ystart zstart xend yend zend
```

In 2D only four values are given:

```
control_print_beam_force_moment_coordinates <index> xstart ystart xend yend
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index`   | Index of the control record. Must match the `control_print_beam_force_moment` index. |
| `xstart ystart [zstart]` | Start point of the cut segment. |
| `xend yend [zend]` | End point of the cut segment. |

The length of the record must be exactly `2 * number_of_space_dimensions`
(4 in 2D, 6 in 3D); anything else, and a zero-length cut, are errors.

## Example

A vertical cut crossing a horizontal beam at `x = 0`:

```
control_print_beam_force_moment 5 -separate_index
control_print_beam_force_moment_coordinates 5 0 -0.5 0 0.5
```

## Differences with the Professional version

- The coordinates record is mandatory for the print to run: without it
  the print reports an error (there is no sensible default cut for
  beams).
