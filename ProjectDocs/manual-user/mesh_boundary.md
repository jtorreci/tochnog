# mesh_boundary + element_boundary + node_boundary

## Description

`mesh_boundary -yes` (manual Professional 6.816) makes Tochnog determine
the boundary of the mesh and fill two print records:

- `node_boundary <node> -yes` (6.877) for every node on the boundary;
- `element_boundary <element> -yes` (6.417) for every element on the
  boundary.

The boundary is the set of element sides (2D edges, 3D faces, 1D end
points) that no other element shares. A node is on the boundary when it
belongs to such a side; an element is on the boundary when it owns such a
side. Interior elements and nodes get no record.

The default is `-no`: without `mesh_boundary -yes` no record is filled
(except `node_boundary` for the element-target contact, see
`contact_target_element_group`). Both records are "for printing only";
the manual says they should not be set by the user.

The records are filled at the start of the first iteration and again
after every mesh change (refinement, element deletion, ...).

## Usage

```
mesh_boundary -yes
```

## Parameters

- `switch`: `-yes` fills the records, `-no` (default) does not.

## Example

A 4 x 4 quad4 mesh (`mesh_boundary1.dat` in `validation-suite/test-2026`):

```
mesh_boundary -yes
```

gives `node_boundary` on the 16 outer nodes and `element_boundary` on the
12 elements with a free side; the 4 interior elements get no record.

## Limitations

Only element types with a side table are scanned: bar2 (1D); tria3,
tria6, quad4, quad9, quad16 (2D); tet4, hex8, hex27 (3D). Other element
types (prism, tet10, hex20, truss, beam, spring, interface, contact) give
no boundary and do not hide one.

## Differences with the Professional version

- Default: the Professional 02-08-2026 fills `node_boundary` and
  `element_boundary` in every run, also without `mesh_boundary`; only
  `mesh_boundary -no` stops it. The GNU follows the manual: default `-no`
  (owner decision 2026-11-02; DIVERGENCES MESH-BOUNDARY-DEFAULT, class A).
  `node_boundary` is also read by the calculation (`distri`, contact,
  `control_new_mesh`, the geometry entity `-node_boundary`), so filling it
  always would change existing models.
- The Professional writes `element_boundary -yes` for EVERY element,
  interior elements included (4 x 4 quad4 mesh: all 16). The GNU marks
  only the elements with a free side, as the manual says ("located on the
  boundary of the mesh") (DIVERGENCES ELEMENT-BOUNDARY-ALL, class D).
- `node_boundary` is the same in both programs.
