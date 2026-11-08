# control_mesh_duplicate_element_group

## Description

Duplicates the elements of a group (Professional manual 6.177): every
element of `element_group_old` gets a copy with the same nodes and the
same element data in group `element_group_new`. Useful, for example, to
put a second material or an interface on the same nodes.

## Usage

```
control_mesh_duplicate_element_group index element_group_old element_group_new
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control index; the copy is made when that index is processed (also without `control_timestep`). |
| `element_group_old` | Group to copy (elements without an `element_group` record are group 0). |
| `element_group_new` | Group of the copies. |

The copies are numbered after the highest element number, in element
order (as the Professional 02-08-2026).

## Example

```
element 2 -quad4 1 2 4 5
element 5 -quad4 2 3 5 6
element 7 -quad4 2 3 5 6
element_group 5 1
element_group 7 2
control_mesh_duplicate_element_group 10 1 4
```

Element 5 gets the copy `element 8 -quad4 2 3 5 6` in group 4. Tests:
`validation-suite/test-2026/duplicate_element_group1.dat`,
`duplicate_element_group2.dat`.

## Differences with the Professional version

None found (the GNU `.dbs` has no `element_group` lines for the default
group 0).
