# node_force / node_inertia / node_slide / node_pressure overrides

## node_force

Discrete nodal force (manual Professional 6.887):
`node_force <index> <fx> <fy> <fz>` (components up to ndim). An
external load applied to the node every step: a positive component
pushes the node towards +x/+y/+z, as `force_point` and
`force_element_edge` and as the Professional. Requires
`materi_velocity`.

Sign fixed 2026-10-06 (BACKLOG C28): before, the GNU applied the force
in the opposite direction (1D bar, `node_force 2 1.0`: velx −2 instead
of +1; the "doubling" was the reversed load compressing the bar to
exactly zero length, a singular state). Verified against the Professional
02-08-2026 on 1D (bar2: velx +1), 2D (quad4: disx 2) and 3D (hex8:
disx 1) probes — identical. Inputs written for the old GNU sign must
flip their `node_force` values.

## node_inertia

Output record (manual 6.889): declare `node_inertia <index> <0.> <0.>`
and tochnog fills it every step with the calculated inertia term
`m*(a + g)` per direction (mass acceleration plus the gravity share).
Typical use — d'alembert: replace the mass inertia by static forces for
the remainder of the calculation with
`control_data_copy index -node_inertia -node_force` +
`control_data_copy_index_factor index -1.` (then remove node_mass).

## node_static_pressure / node_dynamic_pressure / node_total_pressure

User overrides of the calculated pressures (manual 6.886/6.894/6.898):
when the record exists for a node, `post_calcul -static/-dynamic/-total`
uses that value instead of the one derived from the phreatic level.
Requires `groundflow_pressure`.

## node_slide

Additional slide membership (manual 6.893): `node_slide <index> <slide_number>`
makes the node a member of that sliding geometry EVEN IF the node is not
inside the `slide_geometry` itself; the plane normal still comes from
`slide_geometry`. GNU extension note: membership is additive (geometry
OR record), which also covers the Professional use.

## Example

```
node_mass 3  0. 2.0
node_force 3  0. -30.
node_inertia 3  0. 0.
```

Test `node_force_inertia`: with gravity (0,-10) the inertia record of
node 3 fills to ~m*g = -20 (-20.109 with the downward discrete force
active since C28; its target was retargeted on 2026-10-06 from
-19.83 +- 0.2, calibrated with the old reversed sign that pushed the
node up, to -20.11 +- 0.2). Test
`test-2026/node_force_sign1`: unit quad4 bar, `node_force` (1, 0) at
the two free nodes gives disx 2 (= Professional).

## Differences with the Professional version

- `node_inertia`: the Professional record holds one value per dof (16
  in a 2D velocity/displacement/stress model) and rejects `node_mass`
  next to it; the GNU record holds ndim values (m*(a + g)). Test `node_pressure`: `node_static_pressure 1 7.5` makes
`post_calcul -static` return 7.5 for node 1 (calculated: -2.0).
