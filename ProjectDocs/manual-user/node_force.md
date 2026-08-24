# node_force / node_inertia / node_slide / node_pressure overrides

## node_force

Discrete nodal force (manual Professional 6.887):
`node_force <index> <fx> <fy> <fz>` (components up to ndim). Applied
directly to the node rhs each step, same sign convention as
`force_point`/`force_element_edge`. Requires `materi_velocity`.

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
node 3 fills to ~m*g = -20 (exactly -19.83 with the discrete force
active). Test `node_pressure`: `node_static_pressure 1 7.5` makes
`post_calcul -static` return 7.5 for node 1 (calculated: -2.0).
