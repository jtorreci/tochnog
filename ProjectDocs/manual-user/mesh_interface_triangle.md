# mesh_interface_triangle_coordinate / mesh_interface_triangle_element_group / control_mesh_interface_triangle

## Usage (manual Professional 6.856 / 6.857 / 6.201)

Generate zero-thickness interface elements in a 3D mesh with `-tet4`
elements by cutting the mesh with a triangulated plane:

```
mesh_interface_triangle_coordinate  index  x0 y0 z0  x1 y1 z1  x2 y2 z2  [ ... ]
mesh_interface_triangle_element_group index element_group
control_mesh_interface_triangle index -yes
```

- `mesh_interface_triangle_coordinate`: the plane of the interface,
  given as one or more triangles (9 coordinates per triangle). The
  spelling is the SINGULAR `coordinate` — the spelling the Professional
  binary accepts and echoes (the manual text/TOC write
  `mesh_interface_triangle_coordinates`; the corpus file `interface11`
  and the Professional .dbs use the singular form).
- `mesh_interface_triangle_element_group`: the element group attributed
  to the generated interface elements (a `group_interface` group).
- `control_mesh_interface_triangle`: activates the generation at the
  given control index.

Example (corpus test `interface11`): the mesh of three `-tet4` filling
the triangular prism `(0,0,0)-(1,0,0)-(0,1,0)-(0,0,1)-(1,0,1)-(0,1,1)`
is cut at the plane `z = 0.6`:

```
group_type 1  -materi
group_interface 1  -yes
group_interface_materi_memory 1  -total_linear
group_interface_materi_elasti_stiffness 1  1.e11 0.5e11 0.5e11

mesh_interface_triangle_coordinate 10  0. 0. 0.6  100. 0. 0.6  0. 100. 0.6
mesh_interface_triangle_element_group 10  1
control_mesh_interface_triangle 10 -yes
```

## What the generation does

For every `-tet4` crossed by the plane (and whose cut polygon lies
inside one of the triangles of the record):

- the crossing points on the tet edges are created twice (one node per
  side of the interface, same coordinates);
- one zero-thickness interface element is inserted, numbered first:
  triangular cut -> `-prism6` `{side1 x3, side2 x3}`;
  quadrilateral cut -> `-hex8` `{side1 x4, side2 x4}` — with
  `element_group` = the `mesh_interface_triangle_element_group` record;
- the two halves of the tet are re-tessellated and replace the original
  element (1 vertex on a side: `-tet4` + `-prism6`; 2+2: two `-prism6`),
  keeping the `element_group` of the original tet.

For the corpus `interface11` (3 tets -> 9 elements, 18 nodes) the
Professional .dbs numbers the generated elements as interfaces 4,5,6
(tet order) and halves 7..12; the GNU implementation mirrors that
layout.

## Verification

- The generation itself is verified against the Professional .dbs
  structure of `interface11` (element types, numbering order, cut-node
  duplication, element groups).
- **Known blocker (not in this sprint)**: the physics of `interface11`
  does not close yet — the GNU interface element is not consistent with
  `-prism6`/`-tet4` volume neighbours (see the developer manual). The
  corpus test therefore still ends rc=1 on its target.
