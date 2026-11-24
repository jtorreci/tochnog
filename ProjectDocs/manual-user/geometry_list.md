# geometry_list

## Description

`geometry_list` (manual Professional 6.531 and section 9.5) defines a list
of **node numbers** that can be used as a geometrical entity wherever a
geometry is accepted (`bounda_dof`, `bounda_force`, `control_mesh_delete_geometry`,
`data_put`/`data_delete`, `geometry_set`, `control_print_node_geometry`, ...).

A node is "in the geometry" if and only if its number appears in the list.
No coordinates, tolerances or projections are involved, so it is the exact
way to select a known set of nodes (for example nodes exported from a
pre-processor) without building a geometric shape around them.

## Usage

```
geometry_list <index> <node_0> <node_1> ...
... -geometry_list <index> ...
```

| Parameter | Meaning |
|-----------|---------|
| `index` | index of the list (like any geometry index) |
| `node_i` | node numbers (integers), any count up to the record size |

## Semantics

- The list always holds node numbers.
- Element selections (for example `control_mesh_delete_geometry`) use the
  node membership exactly as for any other geometry: an element is selected
  when its nodes are in the list. There is **no** element-number
  interpretation. With a 3-element `-quad4` strip, a list `{3}` deletes
  nothing, a list `{3 4 7 8}` deletes element 3 (its four nodes).
- Evaluations at a coordinate that is not a node (for example integration
  points of contact or area searches) are never "in" a `geometry_list`.
- A node in a list is reported in `node_geometry_present` as
  `-geometry_list <index>` (in alphabetical record order, before
  `-geometry_point`).
- `change_geometry` cannot be applied to a `geometry_list` (the run stops
  with a "Sorry" message); `control_mesh_adjust_geometry` leaves list nodes
  in place (the projection of a node is the node itself).

## Example

```
geometry_list 1 1 5
geometry_list 2 1 2 3 4
geometry_list 6 8

bounda_dof 0 -geometry_list 1 -velx
bounda_dof 1 -geometry_list 2 -vely
bounda_dof 3 -geometry_list 6 -temp
bounda_time 3 0. 0. 1. 10.
```

Self-checking inputs: `validation-suite/test-2026/geometry_list1.dat`
(boundary conditions), `geometry_list_delete1.dat` / `geometry_list_delete2.dat`
(element deletion by node membership). Their content is reproduced in the
developer manual page.

## Differences with the Professional version

None measured. Behaviour verified against the Professional binary
02-08-2026 with the tests above.
