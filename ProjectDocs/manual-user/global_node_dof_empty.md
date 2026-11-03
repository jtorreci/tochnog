# global_node_dof_empty

## Description

Input record (manual Professional 6.552):

```
global_node_dof_empty switch
```

- `-yes` (default): the dofs of a node that is only connected to empty
  elements are set to zero.
- `-no`: they are not set to zero; the node keeps its values.

An element is empty when its group has `group_type` `-empty` (or only
`-empty` / `-none` types), or when its group has a `group_time` record and
is not alive at the end of the time step. The zeroing is done at the end of
every time step. Nodes that belong to no element at all are not touched.

The displacement is not set to zero: it keeps the value of the start of the
step (the zero velocity adds no increment), as in the Professional (see
below).

## Why it is useful

A node that only touches empty (not yet built, excavated, dead) elements has
no stiffness and no material; its stresses, strains, temperatures, ... would
otherwise keep stale or input values in the output and in plots.

## Parameters

- `switch`: `-yes` or `-no`.

## Example

`validation-suite/test-2026/global_node_dof_empty1.dat`: two quad4, the
second one in group 1 with `group_type 1 -empty`; nodes 3 and 6 belong only
to the empty element and have `node_dof` 1.1 .. 2.9 in the input. After the
step:

```
node_dof  3  0 0 1.3 1.4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

(velocity 0, displacement kept, strains / stresses / kappa / work /
temperature 0), equal to the Professional 02-08-2026.
`global_node_dof_empty2.dat` (`-no`) keeps the input values;
`global_node_dof_empty3.dat` uses a dead `group_time` group instead.

## Differences with the Professional version

- The manual says that "the dof's in a node are set to zero"; the
  Professional 02-08-2026 keeps the displacement. The GNU keeps it too
  (project-owner decision 2026-11-03: a zeroed displacement would put the
  nodes of an element that becomes active again back at their start
  position, possibly inside the deformed neighbouring elements;
  DIVERGENCES GLOBAL-NODE-DOF-EMPTY-DISPLACEMENT, class A).
- With `-no` the Professional still gives such a node zero velocity (and
  no displacement increment); the GNU leaves the velocity of a node without
  stiffness unchanged, so the displacement keeps growing with it
  (DIVERGENCES NODE-NO-STIFFNESS-VELOCITY, class B, pre-existing).
- Groups excluded by `mesh_element_group_apply` /
  `control_mesh_element_group_apply` are not empty: they keep their state
  (the Professional ignores that record outside index 0).
- Elements found empty by `materi_density` / `materi_diffusion`
  (`element_empty -yes`) are not counted (the zeroing would remove the
  density of the region that is being filled); not compared with the
  Professional.
