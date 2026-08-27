# control_print_beam_force_moment

## Description

`control_print_beam_force_moment` prints the forces and moments of the
beam / truss / truss-beam elements crossed by a cut segment (manual
Professional 6.262). The cut goes from `(xstart, ystart, zstart)` to
`(xend, yend, zend)` as given by
[`control_print_beam_force_moment_coordinates`](control_print_beam_force_moment_coordinates.md)
(in 2D only `x` and `y` are needed). The forces and moments are written
to the file `beam_force_moment.<index>`.

If the element contains a truss (either a truss element or a truss-beam
element), the truss force is used for the axial force.

The switch must be `-separate_index` (file `beam_force_moment.<index>`)
or `-separate_sequential` (file `beam_force_moment.<seq>`, one file per
print call).

## Input syntax

```
control_print_beam_force_moment <index> <-separate_index | -separate_sequential>
control_print_beam_force_moment_coordinates <index> xstart ystart [zstart] xend yend [zend]
control_print_beam_force_moment_switch <index> -yes | -no
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index`   | Index of the control record. Must match an active `control_timestep`. |
| `switch`  | `-separate_index`: file named with the index. `-separate_sequential`: file named with a running sequential number. |
| coordinates | The cut segment: `xstart ystart zstart xend yend zend` (3D) or `xstart ystart xend yend` (2D). A zero-length cut is an error. |

## Output

One line per crossed element, in the order the cut traverses them (the
lines are sorted by ascending distance from the cut start point):

```
distance force_x_first_node force_y_first_node force_z_first_node
moment_x_first_node moment_y_first_node moment_z_first_node
force_x_second_node force_y_second_node force_z_second_node
moment_x_second_node moment_y_second_node moment_z_second_node
```

The components are given in the LOCAL beam axes: local x is the element
direction (from the first to the second node of the element record),
local y is the in-plane perpendicular of the 2D beam, local z is out of
the beam plane. Because the beam element is a 2D element, the
out-of-plane components (`force_z`, `moment_x`, `moment_y`) are
identically zero.

An element is printed when the minimum distance between its axis and
the cut segment is smaller than `1.e-6 * max(1, cut_length)`.

If no element crosses the cut, **no file is written**. The file is
appended at every print call (one line per element per call).

## Truss elements

- `-truss`: only the axial columns are non-zero: `+N` at the first
  node and `-N` at the second node, where `N` is `ELEMENT_TRUSS_FORCE`
  (positive = tension).
- `-trussbeam`: the axial columns come from the truss force (same
  `+N` / `-N` pattern); the transverse forces and the moments come from
  the beam element.
- `-beam`: no truss; the axial columns are zero (the pure 2D beam has
  no axial stiffness).

## Example

A 2D cantilever beam (1 long, E = I = 1, tip load F = 1e-2) with a cut
crossing the beam at the fixed end:

```
control_print_beam_force_moment 5 -separate_index
control_print_beam_force_moment_coordinates 5 0 -0.5 0 0.5
control_timestep 5 1.e-1 1.e-1
```

Produces `beam_force_moment.5`:

```
0.5 0 -0.00999993055578 0 0 0 -0.00999997222228 0 0.00999993055578 0 0 0 0
```

The first column is the distance from the cut start (0.5). The values
are the element internal nodal forces (ELEMENT_BEAM_MOMENT): transverse
force -F at the fixed node, +F at the tip, bending moment -F*L at the
fixed node, 0 at the tip, axial 0.

## Differences with the Professional version

- The printed components follow the sign convention of the element
  internal force vector (`ELEMENT_BEAM_MOMENT`), verified empirically.
  Use [`control_print_beam_force_moment_switch`](control_print_beam_force_moment_switch.md)
  with `-yes` to invert all signs.
- Selection criterion (documented decision): minimum 3D distance
  between the element axis segment and the cut segment smaller than
  `1.e-6 * max(1, cut_length)`.
- If no element crosses the cut, no file is written (documented
  decision).
- The file is appended at every print call (the manual does not state
  the open mode).
- Values smaller than `1.e-6 * (1 + max|component|)` are snapped to
  zero to remove solver residual noise (documented decision).
- In 3D the local axes follow the beam plane (`group_beam_plane`,
  default `-x -y`), the same frame `beam.cc` uses.
