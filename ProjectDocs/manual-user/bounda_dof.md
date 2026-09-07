# bounda_dof

## Description

Prescribes values (Dirichlet boundary condition) to the degrees of freedom
(dofs) of the nodes selected by `node_range`. It is the professional name of
the mechanism that in the GNU version is called `bounda_unknown`: both
keywords are equivalent and can be used interchangeably.

The prescribed dofs are the primary dofs of the model (`-velx`, `-vely`,
`-temp`, `-pres`, ...). The value each dof takes over time is given with a
matching `bounda_time` record.

For a given boundary index, use either `bounda_dof` (Dirichlet: prescribed
values) or `bounda_force` (Neumann: prescribed forces), but never both.

## Usage

```
bounda_dof <index> <node_range> <dof_0> <dof_1> ...
bounda_time <index> <time_0> <value_0> <time_1> <value_1> ...
```

`node_range` selects which nodes receive the prescribed values:

- a range of node numbers, e.g. `1 2 3` or `-range 1:100`;
- `-all` for every node;
- a geometry, e.g. `-geometry_line 1`, `-geometry_point 1`, `-geometry_surface 1`
  — all nodes lying on that geometry receive the values;
- `-node_set <set>` for a node set.

## Parameters

| Parameter   | Meaning                                                            |
|-------------|--------------------------------------------------------------------|
| `index`     | Boundary index; must match the `index` of the `bounda_time` record.|
| `node_range`| Nodes to apply the boundary to: range, `-all`, `-geometry_*`, `-node_set`. |
| `dof_0 ...` | Primary dofs to prescribe (`-velx`, `-vely`, `-temp`, `-pres`, ...).|

Values over time are set by a `bounda_time <index> <time> <value> ...` record.

## Example

Prescribe zero velocity (`-velx`) at time 0 on the nodes of geometry line 1,
and velocity 1 at time 1:

```
bounda_dof 0  -geometry_line 1 -velx
bounda_time 0  0. 0. 1. 1.
```

The same record can be written with the GNU name:

```
bounda_unknown 0  -geometry_line 1 -velx
bounda_time 0  0. 0. 1. 1.
```

**Note:** `bounda_dof` and `bounda_unknown` are aliases of the same mechanism;
only one of `bounda_dof`/`bounda_unknown` or `bounda_force` may be used per
boundary index.

## Prescribing the total pore pressure: `-topres`

For groundflow models (`groundflow_pressure` initialized) the dof `-pres` is
the hydraulic head solved by the storage equation. If you want to prescribe
directly the TOTAL pore pressure (`p_total`, the pore pressure of the
geotechnical total-stress split) instead of the head, use the dof label
`-topres` (manual Professional 2.4.1; the corpus spelling `-tpres` is
accepted as an alias):

```
bounda_dof 10  -geometry_line 10 -topres
bounda_time 10  -20.
```

The value of `bounda_time` is interpreted as the total pore pressure and
converted per node to the dof value (the DYNAMIC pressure) the groundflow
machinery needs. The conversion depends on the source of the static
pressure at the node (measured against the Professional binary 25-10-2023):

- `post_calcul_static_pressure_height` region covering the node:
  `pres_dof = p_total - p_static` (the head follows the reference height).
- `groundflow_phreatic_level` covering the node:
  `pres_dof = p_total - rho*g*level` (uniform under a constant level; the
  head keeps the load and the reported total pressure follows
  `p_total = h - rho*g*z`, the profile of the manual 2.2.7 undrained
  strategy).
- without any level: `pres_dof = p_total + rho*g*z`.

tochnog sign conventions: g and z negative below the datum.

Example (ground13 of the corpus): bottom nodes `-topres -20` at z=0 and top
nodes `-topres -10` at z=1 with `force_gravity 0. -10.` and
`groundflow_density 1.` give the hydraulic head -20 everywhere (hydrostatic,
no flow) and `-to_pres` = -20/-10 on the bottom/top rows — identical to the
Professional. Example with a phreatic level above the mesh (undrained1 of
the corpus): `-topres` load 0 on `-all` with `groundflow_phreatic_level 1.`
over a mesh at z in [-1,0] gives pres_dof = +10 everywhere (the dynamic
pressure) and `-to_pres` = 0 at z=0 / -10 at z=-1, exactly the
Professional.

## Zero velocity normal to a wall: `-veln`

For velocity dofs, `-veln` prescribes that the nodes do NOT move in the
direction NORMAL to a plane (no-penetration / frictionless-wall condition:
the normal component of the velocity is zero, the tangential components
stay free). Manual Professional 6.22.

```
bounda_dof <index> -geometry_set 1 -veln
```

The normal direction comes from:

- the geometry entity of the record when the nodes are selected through one
  (a `geometry_line` for a wall; for a `geometry_set`, each node uses the
  normal of the FIRST entity of the set that contains it — geometry() first
  match, verified against the Professional binary on validation_8 corner
  nodes lying on two entities);
- a `bounda_normal <index> ...` record when the nodes are selected through a
  node range (the manual 6.22 requires the normal in that case).

`materi_velocity` must be active (`-veln` needs the velocity dofs; validated
at input). The `bounda_time` of the record is irrelevant (manual 6.22).

Internally Tochnog generates multi point constraint records
(`mpc_node_number`/`mpc_node_factor`, marked `mpc_from_bounda -yes`) that
impose the condition:

- one mpc record per boundary node;
- slave dof = the FIRST velocity axis with a non-zero normal component
  (x, y, z order), masters = the remaining axes with a non-zero component,
  factor = `-n_master/n_slave` with `n` the (unit) normal: the record
  `mpc_node_number k <node> -velx <node> -vely` + `mpc_node_factor k <f>`
  imposes `velx = f*vely`, i.e. `n_x*velx + n_y*vely = 0`;
- a zero-component master is omitted: on an axis-aligned wall the record
  degenerates to `mpc_node_number k <node> -vel?` alone (the normal dof
  bounded to zero, e.g. `-vely` on a horizontal wall, `-velx` on a vertical
  one);
- the tie is eliminated inside the system solve (energy-consistent slave
  elimination, the mpc_linear_quadratic mechanism) and the generated
  records are re-created when the mesh changes (refinement, ...).

Record layout measured against the Professional binary (25-10-2023):
validation_8 oblique wall of slope -0.1 generates
`mpc_node_number k <node> -velx <node> -vely` with factor -10; a wall of
slope +0.1 gives factor +10, slope +2 factor +0.5, a vertical wall the
masterless `-velx`, a horizontal wall the masterless `-vely`.

Example (validation_8 of the corpus — continuous extrusion through a
stepped die with frictionless walls):

```
start_define
  lower_edge geometry_line 3
end_define
lower_edge  -1. 0. 2. 0. 1.e-4
...
geometry_set 1  -lower_edge -upper_left_edge -upper_oblique_edge -upper_right_edge

bounda_dof 2 -geometry_set 1 -veln
```

The nodes of the horizontal wall portions get `vely = 0` bounded; the nodes
of the oblique portion `(0,1) -> (1,0.9)` get the tie `velx = -10*vely` (the
flow slides along the wall). A node whose dof is ALSO prescribed by an
explicit `bounda_dof`/`bounda_time` record keeps the direct prescription
(the generated tie is inert there — measured on the Professional binary at
the inflow corner of the wall).

Implementation notes (developer): `bounda.cc::bounda_veln_mpc()` generates
the records (fingerprint bookkeeping per bounda record index, the
mpc_linear_quadratic pattern); `mpc.cc` consumes them (masterless records
bound the slave to zero; the marked records register for the energy
consistent elimination and the direct-bounda slave protection).

## See also
