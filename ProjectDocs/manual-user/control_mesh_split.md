# control_mesh_split

## Description

Splits elements into simpler elements (Professional manual 6.240-6.243):
each `-quad4` into four `-tria3` around a new centre node, each `-hex8`
into twelve `-tet4` around a new centre node (or two `-prism6`), each
`-quad9` into four `-tria6`, each `-hex27` into six `-tet10` and each
`-tria6` into four `-tria3`. Useful to turn a structured mesh into
triangles / tetrahedra or prisms. The split runs in the control index of
the record, also one without `control_timestep`; the mesh is renumbered
afterwards.

## Usage

```
control_mesh_split              index -yes
control_mesh_split_element_from index name
control_mesh_split_element_to   index name
control_mesh_split_only         index geometry_entity geometry_entity_index
```

## Parameters

| Record | Meaning |
|--------|---------|
| `control_mesh_split` (6.240) | `-yes` splits the elements. |
| `control_mesh_split_element_from` (6.241) | Split only the elements with this name (for example `-tria6`); the others stay. |
| `control_mesh_split_element_to` (6.242) | For a `-hex8`: `-tet4` (default) or `-prism6`; for a `-hex27`: `-tet10` (default). |
| `control_mesh_split_only` (6.243) | Split only the elements with all nodes on the geometry. |

The `-tria6` nodes `n0 .. n5` give the `-tria3` `n0 n1 n3`, `n1 n2 n4`,
`n3 n4 n5` and `n1 n3 n4`; the `-hex8` nodes `1 .. 8` give the `-prism6`
`1 2 3 5 6 7` and `2 4 3 6 8 7` (as the Professional 02-08-2026).

## Example

The Professional `test/other/prism6`:

```
number_of_space_dimensions 3
...
element 1  -hex8 1 2 3 4 5 6 7 8
...
control_mesh_split              20  -yes
control_mesh_split_element_to   20  -prism6
control_timestep             30  1. 1.
```

Self-checking inputs in `validation-suite/test-2026`:
`split_element_to_prism6_1`, `split_element_from1`, `split_tria6_1`,
`split_error1`.

## History

F5 T4 U5b (2026-11-08): `_element_from`, `_element_to -prism6` and the
`-tria6` split (before, a `-tria6` was not split). `post_point` and the
other point searches now work in `-prism6` elements (before: an abort,
BACKLOG C125).

## Differences with the Professional version

- `control_mesh_split_element_to -prism18` (hex27 to 18-node prisms) is
  not available: the GNU has no `-prism18` element (the word is unknown).
- A `control_mesh_split_element_to` name other than `-tet4`, `-prism6` or
  `-tet10` is an input error in the GNU; the Professional ignores it
  (DIVERGENCES SPLIT-INPUT-CHECKS).
- The GNU `.dbs` has no `element_group` lines for the default group 0
  (the Professional prints them).
