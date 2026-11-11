# control_mesh_remove_really

Family: `control_mesh_remove_really` (6.232), `_really_activate_all`
(6.233), `_really_activate_factor` (6.234), `control_mesh_remove_keep_
geometry` (6.231). F5 T4 U9e.

## Files and functions

- `mesh_remove.cc` — `mesh_remove(long int icontrol)` (moved from
  `delete.cc`), called from `top.cc` `step_start()` on every step of an
  index with `CONTROL_MESH_REMOVE` (also `step_start(NO)` of an index
  without `control_timestep`).
  - `mesh_remove_dof_scaled()`: the dofs the factor acts on, the stress
    block `stres_indx .. stres_indx+6*nder` and the pore pressure
    `pres_indx .. pres_indx+nder`.
  - `mesh_remove_nodes_in_geometry()`: every node of an inactive element
    in a geometry, `geometry()` with `NODE_START_REFINED` and
    `PROJECT_EXACT` (the test of `element_in_geometry -all`): on the node
    when it exists, else on its kept `node_start_refined` coordinates.
  - `mesh_remove_deactivate()`: keeps the element and deletes it
    (`delete_element`).
  - `mesh_remove_forget()`: deletes the three inactive records.
  - `mesh_remove_activate()`: the reactivation of one element.
- `database.cc` — `CONTROL_MESH_REMOVE_KEEP_GEOMETRY` (INTEGER 2),
  `CONTROL_MESH_REMOVE_REALLY`, `_ACTIVATE_ALL` (INTEGER 1),
  `_ACTIVATE_FACTOR` (DOUBLE 1), all CONTROL with `data_required =
  CONTROL_MESH_REMOVE`; GNU records `ELEMENT_REMOVE_INACTIVE` (INTEGER
  `1+MNOL`, the `ELEMENT` record), `ELEMENT_REMOVE_INACTIVE_GROUP`
  (INTEGER 1, only when the element had an `ELEMENT_GROUP`),
  `ELEMENT_REMOVE_INACTIVE_NODE` (DOUBLE `2*MNOL*MDIM`, `NODE` then
  `NODE_START_REFINED` of its nodes), class ELEMENT, no `data_required`
  (the `ELEMENT` record is absent), no `version_all`.
- `tochnog.h` / `tochnog-mod.h` — the seven enums after
  `ELEMENT_MAP_NODE`.
- `makefile`, `makefile - Windows-Mingw64` — `mesh_remove` object.

## Implementation details

`mesh_remove()`:

1. Records of the index: method and groups, `_geometry`, `_keep_geometry`,
   `_really` (default `-YES`; else `db_error`), `_activate_all` (default
   `-NO`), `_activate_factor` (default 1; outside [0, 1] the Professional
   message and exit). Element range = the highest index of `ELEMENT` and
   `ELEMENT_REMOVE_INACTIVE`; a method other than `-method1` stops as
   before.
2. `node_active[]`: the nodes of the elements active at the start of the
   pass (Professional `node_active_last`).
3. Decision per element, active or inactive, in element order: the
   `-method1` test (group `list[1]`, all nodes are nodes of one active
   element of the groups `list[2..]` not removed earlier in this pass,
   which keeps the sequential result of the old in-place deletion) after
   the `_geometry` restriction; an inactive element not removed is kept
   removed by `_keep_geometry` when all its nodes are in it.
4. Removed: an active element is deleted (`_really -yes`) or deactivated
   (`-no`); an inactive one is forgotten (`-yes`) or stays.
5. Inactive and not removed: `mesh_remove_activate()` in element order,
   with `assigned[]` so that a node is set once.
6. `mesh_has_changed()` after every pass with groups (as before) or when
   anything changed; `nod_nod` then deletes the nodes no active element
   uses.

`mesh_remove_activate()`: stops when `ELEMENT` of that number exists
again (an element numbered after the deactivation). Mean over the element
nodes with `node_active` of `NODE_DOF` and `NODE_DOF_START_REFINED` (0
without such nodes). A node not active gets the mean (scaled dofs times the
factor); when the node no longer exists it is PUT again with `NODE`,
`NODE_START_REFINED` (kept coordinates), `NODE_DOF_START_REFINED` (the
mean) and `NODE_DOF`, as `mesh_map.cc` does. An active node is left alone,
or with `_activate_all -yes` gets the factor times the mean in the scaled
dofs only. The element is created again with `create_element(element,
element, ...)` on its clean index (the solved element data were deleted
with it, so the `DATA_ITEM_SIZE` copy buffer that makes `create_element`
unsafe on a solved element is not involved), plus `ELEMENT_GROUP`. With
`options_element_dof -YES` (`global_element_dof_apply -yes`)
`ELEMENT_DOF_INITIALISED` is set to 1 over the zero `ELEMENT_INTPNT_DOF`
of `create_element`, so the integration point stresses, strains and
history start at zero and the other dofs come from the nodes; measured on
the Professional 02-08-2026, that equals its results (5e-11) for the
factors 0, 0.5 and 1, where the factor does not change the result. With
`-no` the integration points come from the nodes.

## Verification

Professional 02-08-2026 probes (2D quad4 column under self-weight, the
removed element on a coincident container of another group):
default `_really` = `-yes`; inactive element marked
`element_remove_active -yes`, `element_empty -empty`, data cleaned, free
nodes bounded; a step without `control_mesh_remove` changes nothing;
reactivation on a step whose remove does not remove the element; factor
and `_activate_all` of the reactivating index; mean over the element's
active nodes, all dofs, stresses and pressure scaled (temperature not);
elements in order, a node set once, zero without active nodes;
`_keep_geometry` of the current index only, all nodes, followed by the
`_really` of the index. Stage node dofs equal (`control_mesh_remove_
really1`, `_chain1`, all 10 dofs, 1e-9).

## Hardcoded parameters / pending refactorings

- The numbers of inactive elements and their nodes are not renumbered by
  `control_mesh_renumber` or a remesh; an element created on the number of
  an inactive element stops the reactivation with an error.
- `-method3` is not implemented (BACKLOG C126), `-method1` matches node
  numbers.
