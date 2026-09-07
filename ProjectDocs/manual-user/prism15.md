# prism15

## Description

`-prism15` is the 15-node quadratic prism (wedge): the quadratic upgrade
of the 6-node `-prism6`. It is the element the Tochnog Professional uses
for quadratic wedge meshes (corpus `prism15.dat`). The element has:

- 6 corner nodes (base triangle + top triangle),
- 3 mid nodes on the vertical edges,
- 6 mid-edge nodes on the triangular faces (3 base + 3 top).

There are **no** mid nodes on the rectangular side faces (an element
with those is the 18-node prism `-prism18`, not supported).

## Node ordering

The ordering follows the Professional (verified against its `.dbs` of
`prism15.dat`):

1. base triangle corners 1, 2, 3 (at z of the base),
2. top triangle corners 4, 5, 6 (same in-plane order, at the top z),
3. mid nodes of the vertical edges 7, 8, 9 (above corners 1, 2, 3),
4. mid-edge nodes of the base triangle 10, 11, 12 (edges 1-2, 2-3, 3-1),
5. mid-edge nodes of the top triangle 13, 14, 15 (edges 4-5, 5-6, 6-4).

Corpus `prism15.dat` writes exactly this order, e.g. for a prism with
base corners (0,0), (1,0), (0,1) at z=0 and top at z=1:

```
element 1  -prism15 1 2 3  4 5 6  7 8 9  10 11 12  13 14 15
```

The in-plane triangles are 6-node quadratic triangles; the z-direction
is quadratic along the vertical edges (the 15-node wedge is the
serendipity element: no side-face mids).

## Usage

Place the elements in the data part like any other element
(`number_of_space_dimensions 3`):

```
element 1  -prism15 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
```

The element is fully usable in the mixed u-sigma formulation
(`materi_velocity` / `materi_displacement` / `materi_stress`) with
`post_point`, `post_node`, boundary conditions on geometries, etc. —
the same paths as `-hex27`/`-tet10`.

## Notes

- The integration rule is the one measured from the Professional binary:
  3 Gauss points along the prism axis x the 7-point degree-5 triangle
  rule (21 integration points). A uniform stress state is integrated
  exactly: corpus `prism15.dat` returns `sigzz = 1.0` at its
  `post_point` (Professional: 1.0, GNU: 1.000000261, target tolerance
  1e-2).
- Element works in 3D only.
- No `group_integration_points -minimal/-maximal` selection yet: the
  rule is fixed at 21 points (like the other quadratic volume elements).

## Related

- [prism6](prism6.md) — the linear wedge (same base/top ordering).
- [tet10](tet10.md), [hex27](hex27.md) — the other quadratic volume
  elements.
- `element_3d.md` (developer) — volume element architecture.
