# group_time

## Description

`group_time index birth_time death_time` (manual Professional 6.769)
makes the elements of element group `index` active only between two
times. An element is active in a time step that ends at `t` when
`birth_time < t <= death_time` (and at `t = 0` when `birth_time = 0`).
Before its birth and after its death the element adds no stiffness, mass
or load. Useful for staged construction: a strut, a layer or a member
that is added (or removed) during the calculation.

## Parameters

| Parameter | Meaning |
|---|---|
| `index` | element group |
| `birth_time` | time after which the group becomes active |
| `death_time` | time after which the group is inactive again (use a large value, e.g. `1.e20`, for "never") |

## Nodes of inactive elements (since 2026-10-06)

A node whose elements are ALL inactive at the current time (the group is
not yet born, or already dead) is not driven by a `bounda_dof` /
`bounda_unknown` record: the prescribed value is not applied to it, and
the node stays at rest until one of its elements is active. A node that
also belongs to an active element follows its boundary condition as
usual. Nodes that belong to no element at all are not affected by this
rule. This follows the Professional 02-08-2026 (the manual is silent).

## Example

```
number_of_space_dimensions 1
materi_velocity
materi_displacement
end_initia
node 1  0
node 2  1
element 1  -truss 1 2
group_type 0  -truss
group_truss_memory 0 -total_linear
group_truss_area 0  1.0
group_truss_elasti_young 0  1.0
group_time 0  1. 1.e20
bounda_dof   0  1 -velx
bounda_time  0  1.e-1
control_timestep 0  5.e-1 2.
end_data
```

The truss is born at `t = 1`. Node 1 has a prescribed velocity 0.1 from
`t = 0`, but it stays at rest until the truss is born; at `t = 2` both
nodes have `disx = 0.1` (Professional test `other/truss8`).

## Differences with the Professional version

- After the DEATH of all its elements, the Professional also brings the
  free and bounded nodes to rest (velocity 0); the GNU no longer applies
  the `bounda_dof` to such a node, but a node keeps the velocity it had
  (BACKLOG C39). Inputs that kill a group and keep loading its nodes may
  therefore differ.
