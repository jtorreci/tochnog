# condif_convection_edge_normal

## Description

Convection boundary condition on element edges (manual Professional
6.65): heat flux `q = h*(Tenv - T)` applied to the temperature dofs of
the edge nodes, with the tangent stiffness `h` added to the system
matrix. This is the Professional name of the legacy GNU keyword
`condif_convection` (same physics and value layout `h Tenv`; both
names are accepted, the records are independent).

Companions (same index):

- `condif_convection_edge_normal_geometry` — selects the area: either a
  geometry entity (e.g. `-geometry_line 1`; the total edge of an element
  must be inside the geometry) or an explicit node list.
- `condif_convection_edge_normal_element` — restrict to the listed
  elements.
- `condif_convection_edge_normal_element_group` — restrict to elements
  of the listed element groups.
- `condif_convection_edge_normal_element_side` — pairs element/side.
- `condif_convection_edge_normal_node` — only the listed global nodes of
  the edge.
- `condif_convection_edge_normal_element_node` — element + local node
  numbers.

`condif_temperature` must be an unknown; 2D/3D only.

Attention: convection on an INTERIOR edge does not cancel (it applies
from both sides), unlike a prescribed normal flux.

## Usage

```
condif_convection_edge_normal <index> <h> <Tenv>
condif_convection_edge_normal_geometry <index> -<geometry_entity> <index>
```

## Example

```
condif_convection_edge_normal_geometry 0  -geometry_line 1
condif_convection_edge_normal 0  1.0 1.0
```

Slab with k=1, T=0 at the far end (L=1): steady state T(edge) = 0.5
when h=1, Tenv=1 (test `condif_convec`). With `_element` restricting to
an element that does not touch the geometry, nothing is applied
(test `condif_convec_el`).
