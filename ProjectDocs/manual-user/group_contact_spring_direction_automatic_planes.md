# group_contact_spring_direction_automatic_planes

## Description

`group_contact_spring_direction_automatic_planes index switch_x switch_y
switch_z` (manual 6.606, index = element group): which coordinate
components the automatic contact spring direction
(`group_contact_spring_direction_automatic -yes`) may have. Default all
`-yes`.

## Effect on the calculation

- A `-no` switch removes that component of the automatic normal direction,
  and the direction is normalized again. With the centroid direction
  (-0.447, -0.894): `-no -yes -yes` gives (0, -1), `-yes -no -yes` gives
  (-1, 0) (Professional 02-08-2026).
- When nothing is left (for example the direction (0, -1) with switch_y
  `-no`) the calculation stops with
  `Error: please specify a valid contact spring direction for element`.
- The record always takes three switches, also in 2D (the Professional
  rejects fewer).

## Related contact spring behaviour (F5 T5 U6)

- Without `group_contact_spring_plasti_friction` the friction is unlimited
  (manual 6.604, "a very large value"): no frictional slip. Before, GNU used
  f = 0 and lost the tangential force and stiffness.
- In 2D the first tangential direction lies in the plane: t1 = (-n_y, n_x).
  `element_contact_spring_force` index 1 is that tangential force.

## Example

```
group_contact_spring_direction_automatic 1 -yes
group_contact_spring_direction_automatic_planes 1 -no -yes -yes
```

Self-checking inputs in `validation-suite/test-2026`:
`contact_spring_planes1`, `contact_spring_planes2`,
`contact_spring_planes_error1`, `contact_spring_tangential1`,
`contact_spring_inclined1`.
