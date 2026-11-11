# control_mesh_remove_really

## Description

Lets `control_mesh_remove` make the removed elements inactive instead of
deleting them (Professional manual 6.231-6.234). An inactive element takes
no part in the calculation; it becomes active again on a later step whose
`control_mesh_remove` no longer removes it. Its nodes then get start
values from the nodes that stayed active.

- `control_mesh_remove_really index -no`: the elements removed by the
  `control_mesh_remove` of this index become inactive. `-yes` (the
  default, also without the record) deletes them, as before.
- An inactive element is tested again on every step of an index with a
  `control_mesh_remove`: removed again, it stays inactive (or is deleted
  when the `_really` of that index is `-yes`); not removed, it becomes
  active. A step of an index without `control_mesh_remove` leaves it as
  it is.
- `control_mesh_remove_keep_geometry index geometry_item_name
  geometry_item_index`: an inactive element with every node in this
  geometry stays removed at this index, also when the `control_mesh_remove`
  no longer removes it. The `_really` of the index then applies: with the
  default `-yes` the kept element is deleted. Only the record of the
  current index acts, and only on elements that are inactive now (an
  element activated again is not "removed once" any more).
- On reactivation (`control_mesh_remove_really_activate_factor index
  factor`, default 1, from 0 to 1, and `control_mesh_remove_really_
  activate_all index switch`, default `-no`, of the reactivating index):
  - a node of the element that no active element used gets the mean of
    the node dofs of the element's nodes that were active: the stresses
    and the pore pressure times the factor, the other dofs (velocity,
    displacement, temperature, ...) as they are. Without active nodes the
    element's free nodes get zero;
  - a node that was active keeps its dofs; with `_activate_all -yes` its
    stresses and pore pressure are set to the factor times the mean too;
  - when several elements are activated in one step they are handled in
    element order and a node takes its values from the first of them.
- The integration points of the reactivated element start from its nodes
  with the default `global_element_dof_apply -yes`: zero stresses and
  strains (as a new element), the other dofs from the nodes. So with the
  default the factor sets the nodal values the element starts from, which
  the first step then replaces (= Professional 02-08-2026). With
  `global_element_dof_apply -no` the element takes its stresses from the
  nodal values, and the factor acts on the calculation.

## Usage

```
control_mesh_remove                       index -method1 group_0 group_1 ...
control_mesh_remove_really                index switch
control_mesh_remove_keep_geometry         index geometry_item_name geometry_item_index
control_mesh_remove_really_activate_factor index factor
control_mesh_remove_really_activate_all   index switch
```

## Parameters

| Record | Meaning |
|--------|---------|
| `control_mesh_remove_really` (6.232) | `-yes` (default): delete the removed elements. `-no`: make them inactive. Another word is an input error. |
| `control_mesh_remove_keep_geometry` (6.231) | Inactive elements with all nodes in this geometry stay removed at this index. |
| `control_mesh_remove_really_activate_factor` (6.234) | Factor from 0 to 1 (default 1) on the mean stresses and pore pressure given to the nodes of a reactivated element. Outside [0, 1]: `Error, CONTROL_MESH_REMOVE_REALLY_ACTIVATE_FACTOR factor should be between 0 and 1.` |
| `control_mesh_remove_really_activate_all` (6.233) | `-yes`: the already active nodes of a reactivated element get the factor times the mean stresses and pore pressure too. Default `-no`. |

All four need a `control_mesh_remove` of the same index.

The `.dbs` holds, for an inactive element, the GNU records
`element_remove_inactive index name nodes`, `element_remove_inactive_group`
and `element_remove_inactive_node` (the coordinates of its nodes); the
element itself and its nodes that no active element uses are not in the
mesh while it is inactive.

## Example

Element 2 lies on element 3 of group 1 and is removed into the inactive
state at 10; element 3 is deleted at 15; at 20 the remove finds no
container any more, so element 2 is active again and its free nodes 7 8
get the mean dofs of nodes 4 5, the stresses times 0.5:

```
control_mesh_remove                        10 -method1 2 1
control_mesh_remove_really                 10 -no
control_timestep                           10 1. 1.
control_mesh_delete_element                15 3
control_timestep                           15 1. 1.
control_mesh_remove                        20 -method1 7 8
control_mesh_remove_really_activate_factor 20 0.5
control_timestep                           20 1. 1.
```

Self-checking inputs in `validation-suite/test-2026` (values and states =
Professional 02-08-2026 unless a difference below says otherwise):
`control_mesh_remove_really1` (the example, all stages),
`control_mesh_remove_factor1`, `_factor0`, `_factor_default1`,
`_activate_all1` (the nodal values right after the reactivation),
`control_mesh_remove_chain1`, `_chain2` (two stacked elements),
`control_mesh_remove_keep1`, `_keep2`, `_keep3`,
`control_mesh_remove_really_yes1`, `control_mesh_remove_default1`, and the
input errors `control_mesh_remove_factor_error1`,
`control_mesh_remove_really_error1`.

## Notes

- `-method1` matches node numbers (BACKLOG C126): an element is removed
  when all its nodes are nodes of a container element, so the container
  shares them and they stay active. The free nodes of the manual's
  reactivation rule occur when the container is gone (as in the example).

## Differences with the Professional version

- The Professional keeps an inactive element in the mesh
  (`element_remove_active -yes`, `element_empty -empty`) and its unused
  nodes fully bounded; the GNU deletes them and keeps the element in the
  `element_remove_inactive*` records (DIVERGENCES REMOVE-INACTIVE-RECORDS).
- At an index without `control_timestep` the Professional stops with
  `Program error detected for data item: element_remove_active` when it
  removes an element, and leaves a reactivated element inactive (it only
  sets the nodal values); the GNU removes and reactivates there as at a
  timestep (DIVERGENCES REMOVE-NO-TIMESTEP).
- `_activate_all -yes`: the Professional also gives the already active
  nodes the mean velocity and displacement (the displacement of a shared
  node jumps); the GNU follows the manual and changes only their stresses
  and pore pressure (DIVERGENCES REMOVE-ACTIVATE-ALL).
- Pore pressure under gravity: the Professional converts the mean
  pressure with the element middle and the node height, `factor *
  (mean - rho g . x_middle) + rho g . x_node` (factor 0 gives `rho g .
  x_node`, not 0); the GNU gives `factor * mean` as the manual says
  (DIVERGENCES REMOVE-PRESSURE-GRAVITY). In the probes the pressure is
  solved again in the step, so the results agree.
- `control_mesh_remove_really` with a word other than `-yes` / `-no` is
  an input error in the GNU; the Professional takes it as `-no`
  (DIVERGENCES REMOVE-REALLY-WORD).
