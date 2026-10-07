# contact_target_geometry / contact_target_geometry_switch

## What it does

`contact_target_geometry index geometry_entity_item geometry_entity_index`
(manual Professional 6.106) defines a contact target geometry. Contacting
nodes are forced to stay on the outward-normal side of that geometry; a
penetrating node gets a normal penalty force
`contact_penalty_velocity * penetration` (plus friction, see
`contact_plasti_friction`).

`contact_target_geometry_switch index -yes` (6.107) flips the sign of the
outward normal of the geometry with the same index.

These are the Professional names of the GNU records `contact_geometry` and
`contact_geometry_switch`. Both spellings are accepted and mean the same.

## Why it is useful

Professional input files (for example the validation case `validation_10`
and the tests `contac2` and `pressed_rubber_cylinder`) use the
`contact_target_*` names, often on an index other than 0. They now run
unchanged in the GNU and give the Professional result.

## Syntax

```
contact_target_geometry        index  -geometry_xxx  geometry_index
contact_target_geometry_switch index  -yes|-no
```

| Parameter | Meaning |
|-----------|---------|
| `index` | Number of the contact geometry; any index >= 0, several may be used. |
| `-geometry_xxx geometry_index` | The geometry entity (point in 1D, line/circle/... in 2D, ...). |
| switch | `-yes` flips the outward normal of the geometry with the same index. Default `-no`. |

If both `contact_geometry` and `contact_target_geometry` are given for the
same index, `contact_geometry` is used.

## Example

1D bar whose free end node 2 moves towards a point at x = 1
(Professional test `other/contac2`, `validation-suite/test-2026/contact_target_geometry_index1.dat`):

```
geometry_point 1  1.3 0.3
node_dof 2  0.03  0. 0. 0. 0. 0. 0.
contact_penalty_velocity  10.
contact_target_geometry 1  -geometry_point 1
```

After one step node 2 stops at x = 1.0 (Professional 1.0).

## Status

Implemented (2026-10-07, K05 / BACKLOG C5). Before, only index 0 was read
and, when it was, the GNU aborted with "Program error detected for data:
contact_geometry. Data is allocated in a parallel loop"; other indices were
silently ignored.
