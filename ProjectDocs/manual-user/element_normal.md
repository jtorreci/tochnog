# element_normal

## Description

`element_normal index normal_x normal_y normal_z` (manual Professional
6.445) is a generated record: "the normal of an element in the inward normal
direction of the mesh boundary. The side of an element is on the mesh
boundary if no other element is connected to that side." In 1D only
`normal_x` is filled, in 2D `normal_x normal_y`, in 3D all three.

It is used by the plane cut-off of the direct plasticity laws when the
plane normal is taken automatically:

- `group_materi_plasti_mohr_coul_direct_normal_automatic` (6.728)
- `group_materi_plasti_tension_direct_normal_automatic` (6.740)

Each law uses its own switch: the tension law and the Mohr-Coulomb law
may have different planes (one explicit, one automatic), see
`group_materi_plasti_mohr_coul_direct.md`.

Typical use: a layer of isoparametric elements along a boundary (a soil
interface layer, a crack band) that should only carry a limited normal or
shear traction on the plane parallel to that boundary.

## How the normal is computed

The rule was measured on the Professional binary (02-08-2026) with about 30 probe
meshes (quad4, quad9, tria3, tria6, tet4, hex8, 1D bar2, distorted and
rectangular meshes, several element groups, trusses, element deletion,
large rotation):

1. A side of an element is a **boundary side** when no other isoparametric
   element holds all its nodes. Element groups do not matter (two groups
   sharing a side are connected); trusses, beams and springs do not close
   a side.
2. The **node normal** is the normalized sum, over the boundary sides of
   the node, of the side normal weighted by the side shape function of the
   node (`integral N_node n dA`). For straight sides this is the side
   normal times the side length (area in 3D) divided by the number of side
   nodes; each side normal points into its element. A node without
   boundary sides has a zero normal.
3. The **element normal** is the normalized sum of the node normals of the
   element nodes, or zero when that sum vanishes (an interior element, or
   an element whose boundary sides cancel, such as the middle element of a
   single row).

Examples (all identical in the GNU and the Professional):

| Mesh | Element | element_normal |
|------|---------|----------------|
| column of 2 quad4 (1 x 1) | lower / upper | (0, 1) / (0, -1) |
| 3x3 quad4 of 2 x 1 cells | corner (bottom left) | (0.60706, 0.79465) |
| same | edge (bottom middle) / interior | (0, 1) / (0, 0) |
| single row of 3 quad4 | middle | (0, 0) |
| column of 2 hex8 | lower / upper | (0, 0, 1) / (0, 0, -1) |

The record is generated when a group needs it (one of the two
`_normal_automatic` records is present, whatever its switch value) and it
is regenerated after every mesh change (for example
`control_mesh_delete_element`), from the current node coordinates. It is
**not** updated for deformation: a body rotated 90 degrees with
`mesh -follow_material` keeps the initial normal until the mesh changes.
The record is printed in the database. A value given in the input file is
discarded and recomputed (as in the Professional).

Only a `_normal_automatic` record switched `-yes` makes the direct
cut-off use the element normal. With the switch `-no` the record is still
generated (and printed), but the cut-off ignores it: the group keeps the
full principal-stress cut-off, exactly as without the `_normal_automatic`
record (as in the Professional; test
`validation-suite/test-2026/element_normal2d_column_switch_no1.dat`, a
column of two quad4 stretched in x: `sigxx` 0 in the upper element with
`-no`, 1 with `-yes` in `element_normal2d_column_switch_yes1.dat`).

## Usage

The user does not set `element_normal`; it appears in the `.dbs` file:

```
element_normal  1  0 1
element_normal  2  0 -1
```

It can be checked with a target (value number 0 = `normal_x`):

```
target_item  0  -element_normal 2 1
target_value 0  -1. 1.e-10
```

## Example

```
element 1  -quad4 1 2 3 4
element 2  -quad4 3 4 5 6
element_group 1 1
element_group 2 2
group_type 2  -materi
group_materi_elasti_young 2  1.0
group_materi_plasti_tension_direct 2  0.0
group_materi_plasti_tension_direct_normal_automatic 2  -yes
```

The upper element gets `element_normal` (0, -1): its zero tension cut-off
acts on the horizontal plane, so pulling the top up opens the element
freely (Professional test `tension_direct1`: `disy` 0 in the middle nodes,
all stresses 0).

Self-checking inputs in `validation-suite/test-2026/`:
`element_normal2d_column1.dat`, `element_normal2d_corner1.dat`,
`element_normal2d_mesh_change1.dat`, `element_normal3d_column1.dat`,
`element_normal3d_distorted1.dat`, `element_normal2d_row_explicit1.dat`,
`element_normal2d_row_layer1.dat`, `element_normal2d_row_mixed1.dat`. All
pass on the Professional 02-08-2026.

## Zero element normal: input error

A zero `element_normal` cannot define a cut-off plane. When an element of
a group with `group_materi_plasti_tension_direct_normal_automatic -yes` or
`group_materi_plasti_mohr_coul_direct_normal_automatic -yes` (and no
explicit `_direct_normal` of that SAME law) gets a zero normal, the
calculation stops with an input error, for example (single row of three
quad4, `validation-suite/test-2026/element_normal2d_row_error1.dat`):

```
Error: element_normal is zero for element : 2
Error: element_group of that element is : 0
The element is not on the mesh boundary (interior element,
or its opposite boundary sides cancel), so its normal cannot
define the direct cut-off plane of : group_materi_plasti_tension_direct_normal_automatic
Restrict the element group to the boundary interface layer,
or give the normal explicitly with : group_materi_plasti_tension_direct_normal
If the element rests on a supported face (e.g. a single row
of elements on a support), give that face normal with : group_materi_plasti_tension_direct_normal
Error detected for data item : group_materi_plasti_tension_direct_normal_automatic
Error detected for record : 0
```

The message names the explicit record matching the automatic one that
triggered it (`_tension_direct_normal` or `_mohr_coul_direct_normal`).
The check is made per law: an explicit normal of the OTHER law does not
define this law's plane, so it does not suppress the error. For example
`group_materi_plasti_mohr_coul_direct_normal 0  0. 1.` together with
`group_materi_plasti_tension_direct_normal_automatic 0  -yes` still stops
on a zero element normal, naming the tension records
(`element_normal2d_row_mixed_error1.dat`); the swapped case on elements
with a nonzero normal runs (`element_normal2d_row_mixed1.dat`).
Typical cases: a group that covers the whole mesh (interior elements), a
single element (all sides are boundary sides and cancel), the middle
elements of a single row. Remedies:

- restrict the `_normal_automatic` group to the boundary layer, e.g. only
  the end elements of a row (`element_normal2d_row_layer1.dat`);
- or give the normal explicitly, e.g. the normal of the supported face of
  a row on a support (`element_normal2d_row_explicit1.dat`):
  `group_materi_plasti_tension_direct_normal 0  0. 1.`

The check runs whenever the record is (re)generated, so it also fires
after a mesh change (for example element deletion) that leaves such an
element with a zero normal. Groups with the switch `-no`, elements of
other groups and non-isoparametric elements (trusses, beams, springs) are
not checked.

## Differences with the Professional version

- An element with a zero `element_normal` that uses `_normal_automatic
  -yes`: both versions stop. The Professional aborts on an internal
  assertion ("Program error, assertion failed ... !array_null( directions,
  MDIM*MDIM )"); the GNU gives the readable input error above (element,
  group, remedies). Until 2026-10-01 the GNU silently applied the full
  principal-stress cut-off instead. Until 2026-10-03 an explicit normal of
  the other law also suppressed the error (the law then silently took the
  principal cap); the Professional also aborts in that mixed case.
- `_normal_automatic -no`: both versions take the full principal-stress
  cut-off (no plane). Until 2026-10-02 the GNU still used the element
  normal with `-no` (plane cut-off).
- Quadratic elements with **curved** sides (quad9/tria6 whose mid-side
  nodes are off the chord): the Professional gives values that do not
  follow any consistent rule found (some node normals are axis-aligned),
  the GNU keeps the consistent rule. Straight-sided quadratic elements
  agree exactly.
- Supported element types for the boundary sides: bar2 (1D), tria3, tria6,
  quad4, quad9, quad16 (2D), tet4, hex8, hex27 (3D). Other types (tet10,
  hex20, prism, ...) contribute no boundary side in the GNU.
- The Professional also prints `node_normal`, `node_boundary` and
  `element_boundary` records in the same situation; the GNU does not
  (the node normals are internal).
- `global_normal -yes` (manual 9.1, generates the normals without a
  `_normal_automatic` group) is not accepted by the GNU.
