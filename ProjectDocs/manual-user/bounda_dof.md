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
