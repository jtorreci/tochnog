# control_mesh_extrude_element

## Description

Manual Professional 6.181: "If you extrude -tria6 elements, you can set
name either to -prism12 or -prism18. Then either the 12 node or 18 node
prismatic elements will be generated. Default, if this
control_mesh_extrude_element is not set, then -prism18 is used for name."

It goes with [control_mesh_extrude](control_mesh_extrude.md): every
`-tria6` of the 2D mesh becomes, per element layer, a
[-prism18](prism18.md) with the 6 nodes of its base plane, of the mid
plane of the layer and of its top plane, each in the tria6 order. The
other 2D elements are extruded as before (`-bar2` to `-quad4`, `-tria3` to
`-prism6`, `-quad4` to `-hex8`, `-quad9` to `-hex27`); the record has no
effect on them.

## Usage

```
control_mesh_extrude          index z0 z1 ...
control_mesh_extrude_element  index -prism18
```

The record needs `control_mesh_extrude` with the same index.

| name | meaning |
|---|---|
| `-prism18` | the default: 18-node prisms |
| `-prism12` | the 12-node prism; **not available yet** in the GNU version (pending, F5 T4 U9c): input error `Error: control_mesh_extrude_element -prism12 is not available yet (the 12-node prism is pending); use -prism18.` |
| anything else | input error `Error: control_mesh_extrude_element must be -prism18 or -prism12.` |

## Numbering

As for the `-quad9`: with a `-tria6` in the mesh every element layer also
has a mid plane, numbered as a plane of its own (plane `p` of node `i` is
node `p*n + i`, layer `l` of element `e` is element `l*m + e`). Nodes of
the mid planes that no element uses (the nodes of linear elements) are
deleted.

Example (two tria6 on the unit square, `control_mesh_extrude 10 0. 1. 3.`,
`control_mesh_extrude_n 10 2 1 1`, so 3 layers and 7 planes of 9 nodes):

```
element 1 -tria6 1 2 3 4 5 6
element 2 -tria6 3 7 8 5 9 6
```

becomes

```
element 1 -prism18 1 2 3 4 5 6 10 11 12 13 14 15 19 20 21 22 23 24
element 2 -prism18 3 7 8 5 9 6 12 16 17 14 18 15 21 25 26 23 27 24
element 3 -prism18 19 20 21 22 23 24 28 29 30 31 32 33 37 38 39 40 41 42
...
```

the same mesh as the Professional 02-08-2026.

## Differences with the Professional version

- A name other than `-prism12` / `-prism18` is an input error in the GNU;
  the Professional takes it as `-prism18` (DIVERGENCES
  EXTRUDE-ELEMENT-NAME, class C).
- `-prism12` is refused by the GNU until the 12-node prism exists; the
  Professional generates it (DIVERGENCES EXTRUDE-ELEMENT-PRISM12, class B).
- A mesh with `-tria6` and a linear element (`-quad4`, `-tria3`, `-bar2`)
  without a `-quad9`: the GNU generates the `-prism18` of the manual (with
  mid planes); the Professional generates `-prism12` elements also with
  `-prism18` (and then stops in the calculation with `Error detected for
  data item : element_intpnt_iso_coord`). DIVERGENCES EXTRUDE-ELEMENT-MIXED,
  class D.

## Tests

`validation-suite/test-2026/`: `control_mesh_extrude_element_default1`
(mesh = Professional), `_prism18_1` (mixed mesh with a quad4, new groups,
self-weight reaction -4 = analytic), `_quad4_1` (no effect on quad4 /
tria3, mesh = Professional), `_error1` (bad name), `_error2` (`-prism12`).
