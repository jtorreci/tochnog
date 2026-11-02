# element_beam_direction

## Description

`element_beam_direction <element> x_x x_y x_z y_x y_y y_z z_x z_y z_z`
is an output record (manual Professional 6.414): the local beam axes of
a beam element, as three unit vectors in global axes, local x first,
then local y and local z. Local x runs from the first to the second
node of the beam; local z comes from `element_beam_direction_z`,
`group_beam_direction_z`, `group_beam_direction_z_reference_point` or the
default (see [group_beam_direction_z](group_beam_direction_z.md));
local y = z x x.

For the classic 2D beam the record is `c s 0 -s c 0 0 0 1` with (c, s)
the beam direction in the x-y plane.

## Example

`beam_3d_default_axes1.dat` (validation-suite/test-2026), a beam from
(0,15,0) to (1,16,1) without a local z record:

```
element_beam_direction 6 0.5774 0.5774 0.5774 -0.7071 0.7071 0 -0.4082 -0.4082 0.8165
```

## Differences with the Professional version

- 2D beam: no Professional reference (the Professional 02-08-2026 2D
  beam stops with a singular matrix on the probes); the GNU writes the
  9 values by analogy with 3D (local z = global z).
- Until 2026-11-02 the GNU wrote 2 values (the in-plane direction of
  the planar beam), BEAM-DIRECTION-LENGTH, resolved.
