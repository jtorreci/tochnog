# control_mesh_split

## Files and functions

- `split.cc` — `mesh_split(version)`: called from `top.cc` (`step_start`
  block of `CONTROL_MESH_SPLIT`, then `renumbering()`).
- `database.cc` — `CONTROL_MESH_SPLIT` (INTEGER switch),
  `CONTROL_MESH_SPLIT_ONLY` (geometry), `CONTROL_MESH_SPLIT_ELEMENT_FROM`
  and `CONTROL_MESH_SPLIT_ELEMENT_TO` (INTEGER, 1, an element name word;
  both require `CONTROL_MESH_SPLIT`; F5 T4 U5b).
- `point_el.cc` — `point_el()`: the `-PRISM15` Newton branch also serves
  `-PRISM6` (`prism15_shape()` evaluates the linear prism for nnol 6, and
  the final weights are `Lk (1 -/+ zeta)/2`).

## Implementation details

`mesh_split()` copies the mesh into `VERSION_TMP` element by element:
an element that is in the `_only` geometry, has the `_element_from` name
(when given) and is a `-quad4`, `-quad9`, `-hex8`, `-hex27` or `-tria6`
is replaced by its split pattern (`split_nodes[]`, `nsplit_element` x
`nsplit_node`, plus `nnew_node` new nodes averaged from `new_node[][]`);
`split_center` marks the patterns whose last node is the new centre node
(quad4 -> tria3, hex8 -> tet4). F5 T4 U5b patterns: `-tria6` -> 4
`-tria3` (0 1 3, 1 2 4, 3 4 5, 1 3 4) and, with `_element_to -prism6`,
`-hex8` -> 2 `-prism6` (0 1 2 4 5 6, 1 3 2 5 7 6), no new nodes. Both
equal the Professional 02-08-2026 (`other/prism6`, probes). The
`_element_to` value is checked once (`-tet4`, `-prism6`, `-tet10`).

Verified: the Professional `other/prism6` (mesh and `post_point` sigzz
1), probes for `_from` on a tria6 + quad4 mesh, the tria6 split and the
invalid name. Differential against the U5a commit (33 inputs with
`control_mesh_split`, `-prism6` or `-tria6`): all `.dbs` identical except
`prism6`, which now runs (before: unknown record).

## Pending

- `-hex27` -> `-prism18` (no prism18 element).
- `force_element_edge` and the other edge / face loads of `area.cc` have
  no `-prism6` side table (mixed triangle / quadrilateral faces): the
  Professional `other/extrude2` stops there (BACKLOG C125).
