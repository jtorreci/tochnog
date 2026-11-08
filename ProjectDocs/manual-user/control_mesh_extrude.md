# control_mesh_extrude

## Description

Extrudes a flat 2D mesh (the x-y plane of the input) to 3D, as the
Professional manual 6.179-6.185 describes. Each 2D element is swept
through a stack of element layers: `-bar2` becomes `-quad4`, `-tria3`
becomes `-prism6`, `-quad4` becomes `-hex8` and `-quad9` becomes `-hex27`.
Useful to build a solid volume from a 2D profile without remeshing.

The model must run in 3D (`number_of_space_dimensions 3`). The extrusion
runs once, when its own control index is processed, also when that index
has no `control_timestep` (the mesh can be built in a control index of
its own).

## Usage

```
control_mesh_extrude                                index z0 z1 z2 ...
control_mesh_extrude_n                              index n0 n1 n2 ...
control_mesh_extrude_direction                      index -y
control_mesh_extrude_element_group_new              index old_0 old_1 ... new_00 new_01 ... new_10 ...
control_mesh_extrude_contact_spring_element_group   index group_0 group_1 ...
control_mesh_extrude_contact_spring_element_group_new index new_0 new_1 ...
```

## Parameters

| Record | Meaning |
|--------|---------|
| `control_mesh_extrude` (6.179) | Coordinates `z0 z1 ...` of the layer planes. The 2D nodes become the plane `z0`; the mesh is extruded up to the last value. |
| `control_mesh_extrude_n` (6.185) | Number of element layers in each interval: `n0` between `z0` and `z1`, `n1` between `z1` and `z2`, and so on. The last value is a dummy (use 1). Default: one layer per interval. |
| `control_mesh_extrude_direction` (6.180) | `-y`: extrude along the global y-axis; the 2D y coordinate becomes z. Default (or `-z`): along z. |
| `control_mesh_extrude_element_group_new` (6.184) | The old element groups of the 2D elements, then for each old group one new group per interval `z(i)..z(i+1)`. A negative new group drops the elements of that interval. Elements of other groups keep their group. |
| `control_mesh_extrude_contact_spring_element_group_new` (6.183) | A `-contact_spring2` is generated between each node of the plane `z0` and its copy on the last plane. With one value, every spring gets that group. |
| `control_mesh_extrude_contact_spring_element_group` (6.182) | The old groups that go with the new groups of `_new`: a spring gets the new group of the first listed old group attached to its node; nodes attached to none of them get no spring. |

## Numbering

As in the Professional 02-08-2026: plane `p` (0 = `z0`) of node `i` is node
`p*n + i`, and layer `l` (0 = the lowest) of element `e` is element
`l*m + e`, with `n` / `m` the highest node / element number (plus 1 when
number 0 is used). Layer 0 keeps the 2D element numbers. Nodes that no
element uses any more (a negative new group) are deleted. The contact
springs are numbered after the extruded elements.

Node order of the new elements (Professional): `-quad4 a b c d` becomes
`-hex8 a b a' b' c d c' d'`, with `'` the copy on the next plane, so an
interface `-quad4` keeps its two sides; `-tria3 a b c` becomes `-prism6
a b c a' b' c'`; `-bar2 a b` becomes `-quad4 a b a' b'`. A `-hex27` has
the 9 nodes of its base plane, of its mid plane and of its top plane in
the quad9 order. With a `-quad9` in the mesh every layer also has a mid
plane, numbered as a plane of its own.

## Example

The Professional `test/other/extrude1`: two old groups, the elements of
group 101 only up to z = 1.

```
number_of_space_dimensions 3
end_initia
node 1  0.0 0.0 0.0
...
element_group -ra 1 2 3 -ra  101
element_group 4  102
control_mesh_extrude                   10  -1. 0. +1. +2.
control_mesh_extrude_n                 10  5 5 5 1
control_mesh_extrude_element_group_new 10  101 102
                                           1 1 -1
                                           2 3 4
end_data
```

15 element layers; group 101 gets group 1 in the first two intervals and
is dropped in the third, group 102 gets 2, 3 and 4: 45 elements and 119
nodes, the same mesh as the Professional.

Self-checking inputs in `validation-suite/test-2026`: `extrude_layers1`,
`extrude_direction1`, `extrude_element_group_new1`,
`extrude_contact_spring1`, `extrude_contact_spring2`, `extrude_single1`,
`extrude_error1`, `extrude_error2`, `extrude_error4`,
`extrude_quad9_shared_nodes1`, `post_point_extrude_hex27_1`.

## History

- F5 T4 U5a (2026-11-08): the layers follow manual 6.179 / 6.185. Before,
  the linear elements kept the 2D nodes as an extra bottom layer at z = 0
  (`control_mesh_extrude 10 0. 1.` gave a collapsed first layer and a
  second one), `control_mesh_extrude_n` was ignored for them, a `-bar2`
  stayed in the mesh next to its `-quad4`, the hex8 node order was
  `a b c d a' b' c' d'`, the numbering was the GNU's own and the
  extrusion ran at the first timestep only. The Professional
  `interface_bar2_hex8` now gives the same mesh and nodal results
  (7.6e-12) as the Professional (BACKLOG C50).
- F4 T12 (2026-10-11, BACKLOG C47): the hex27 inner planes hold one node
  per source node, shared by the neighbouring elements.

## Differences with the Professional version

- A single coordinate (`control_mesh_extrude 0 1.`) extrudes from the 2D
  plane z = 0 to it (GNU extension, the GNU usage before U5a); the
  Professional ignores the record. `control_mesh_extrude_direction` other
  than `-y` / `-z` is an input error in the GNU; the Professional accepts
  it (DIVERGENCES EXTRUDE-INPUT-CHECKS).
- The `.dbs` prints the GNU element name `-contactspring` for the
  Professional `-contact_spring2` (same element).
- `control_mesh_extrude_element` (6.181, `-prism12` / `-prism18` for
  `-tria6`) is not available: the GNU has no 12 / 18 node prism (F5 T4
  U9).
- hex27 nodal stresses and `post_point` in a hex27 extruded mesh: equal
  to the Professional since F4 T6f C70 / C57 (2026-10-13).
