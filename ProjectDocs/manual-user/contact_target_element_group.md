# contact_target_element_group (element-target contact)

## What it does

`contact_target_element_group group_0 group_1 ...` (manual Professional
6.105) makes the elements of the listed element groups act as contact
targets. A node of another element that penetrates a target element is
pushed back through the nearest boundary side of that element with the
penalty `contact_penalty_velocity * penetration`; the reaction goes to
the side nodes, so momentum (and, with the other penalties, pressure or
temperature) is exchanged between contacter and target.

A target side is used only when all its nodes are on the boundary of the
mesh (`node_boundary`, manual 6.875). Since 2026-10-07 the GNU generates
`node_boundary` itself when `contact_target_element_group` is given and
the input has no `node_boundary` record: a node is on the boundary when
it belongs to an element side that no other element shares. The generated
records are printed in the dbs file, as the Professional does.

## Why it is useful

Professional inputs (`contac1`, `contac4`, `axisym3`) rely on the
automatic boundary; before, the GNU found no boundary side and no contact
happened at all (the contacter flew through the target).

## Syntax

```
contact_target_element_group  group_0 group_1 ...
contact_penalty_velocity      penalty
```

| Parameter | Meaning |
|-----------|---------|
| `group_i` | Element groups whose elements are targets (`-bar2`, `-tria3`, `-quad4`, `-tet4`, `-hex8`). |

Hand-written `node_boundary index -yes` records (GNU 2014 inputs such as
sfnet `contac1`) are still accepted; then no boundary is generated.

## Example

1D, Professional test `other/contac1`
(`validation-suite/test-2026/contact_element_target_boundary1.dat`):

```
element 1  -bar2 1 2
element 2  -bar2 3 4
element 3  -bar2 4 5
element_group 1 1
element_group 2 2
element_group 3 2
node_dof 2 +0.03
contact_target_element_group  2
contact_penalty_velocity  1.e10
```

After one step node 2 has velx 0.02 and node 3 velx 0.01 (Professional
0.02 / 0.01). A 2D quad4 case is
`validation-suite/test-2026/contact_element_target_boundary2d1.dat`.

## Differences with the Professional version

- The Professional stops with "a target element has all nodes on the
  boundary of the mesh" when the target is a single layer of elements
  (manual: "more than 1 layer"); the GNU does not check it (BACKLOG C40).
- The Professional writes `node_boundary` and `element_boundary` in every
  run; the GNU only generates `node_boundary`, and only for element-target
  contact.
- hex8 target sides are not detected in the GNU (BACKLOG C40).
