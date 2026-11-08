# control_mesh_renumber_element_group_offset / control_mesh_renumber_element_geometry_offset

## Description

Offsets the element group and element geometry numbers while
renumbering (Professional manual 6.236, 6.237): with
`control_mesh_renumber` at the same index, every element gets its
`element_group` and / or `element_geometry` plus the offset. Useful to
merge meshes whose group numbers would collide.

## Usage

```
control_mesh_renumber                        index lowest_element lowest_node
control_mesh_renumber_element_group_offset    index offset
control_mesh_renumber_element_geometry_offset index offset
```

An element without the record counts as 0, so it ends with the offset
itself (as the Professional 02-08-2026).

## Example

```
element_group 5 1
element_group 7 2
element_geometry 5 3
control_mesh_renumber 10 1 1
control_mesh_renumber_element_group_offset 10 100
control_mesh_renumber_element_geometry_offset 10 7
```

The three elements become 1, 2, 3 with groups 100, 101, 102 and
element_geometry 7, 10, 7. Test:
`validation-suite/test-2026/renumber_offset1.dat`.

## Differences with the Professional version

None found.
