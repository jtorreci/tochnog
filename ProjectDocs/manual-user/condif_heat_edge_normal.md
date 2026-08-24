# condif_heat_edge_normal

## Description

Distributed prescribed heat flux normal to the edge of an element
(manual Professional 6.72). The distributed heat is translated into
equivalent nodal heat on the temperature degrees of freedom of the edge
nodes. Also specify `condif_heat_edge_normal_geometry` to select the
edge; optionally `condif_heat_edge_normal_time` (or `_sine`) for a
temporal load.

Restriction variants (same index):

- `condif_heat_edge_normal_element` — restrict to the listed elements.
- `condif_heat_edge_normal_element_group` — restrict to elements of the
  listed element groups.
- `condif_heat_edge_normal_element_side` — pairs element/side.
- `condif_heat_edge_normal_node` — apply only on the listed global nodes.
- `condif_heat_edge_normal_element_node` — element + local node numbers.
- `condif_heat_edge_normal_element_node_factor` — per-local-node factors
  (first value = element number).
- `condif_heat_edge_normal_factor` — spatial polynomial
  `a0 + a1*x + a2*x^2 + ...` multiplying the heat.

`condif_temperature` must be an unknown.

Attention (manual): only for linear and quadratic isoparametric
elements. If used INSIDE a FE mesh the elements on each side get the
distributed heat with opposite normals, so the total flux normally
cancels.

## Usage

```
condif_heat_edge_normal <index> <heat>
condif_heat_edge_normal_geometry <index> -<geometry_entity> <geometry_index>
```

## Example

```
geometry_line 1  0. 0. 1. 0. 0.

condif_heat_edge_normal_geometry 0  -geometry_line 1
condif_heat_edge_normal 0  0.1
```

A column with conductivity 0.1 and T=0 prescribed at the top receives
q=0.1 through the bottom edge: steady state gives T(bottom) = q*L/k =
2.0 and a linear profile (test `condif_heat_edge`).
