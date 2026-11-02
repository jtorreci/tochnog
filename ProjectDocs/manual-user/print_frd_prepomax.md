# print_frd_prepomax

## Description

`print_frd_prepomax switch` (manual Professional 6.980) is the default of
`control_print_frd_prepomax index switch` (6.290) for every
`control_print_frd` file; a control record overrules it for its index.
Default `-no`. With `-yes` a 2D mesh is extruded to a thickness of 1
(prepomax cannot plot 2D elements, manual 6.288): every node n gets a copy
at z = 1 numbered n + (highest node number + 1) with the same results, a
tria3 / tria6 becomes a penta6 (frd type 2) and a quad4 / quad9 a hex8 (frd
type 1) on the corner nodes. A 3D mesh is unchanged.

## Usage

```
control_print_frd 0 -yes
print_frd_prepomax -yes
```

| Parameter | Meaning |
|---|---|
| switch | `-yes`: extrude 2D meshes for prepomax. `-no` (default). |

Tests `validation-suite/test-2026/print_frd_prepomax1/2`,
`control_print_frd_prepomax1` (quad4: nodes 1 2 4 3 6 7 9 8; tria3: 1 2 3
6 7 8 = Professional 02-08-2026).

## Differences with the Professional version

The Professional writes 0 as element count in the element block header of
an extruded quad9 mesh; the GNU writes the real count. The tria6 extrusion
(penta6 on the corners) is not measured on the Professional. Other frd
layout differences: BACKLOG C27.
