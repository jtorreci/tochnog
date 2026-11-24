# condif_radiation_edge_normal

## Description

Radiation boundary condition on element edges (manual Professional
6.92): heat flux `q = alpha_r*(Tr^4 - T^4)` with tangent stiffness
`4*alpha_r*T^3`, applied to the temperature dofs of the edge nodes.
This is the Professional name of the legacy GNU keyword
`condif_radiation` (same physics and value layout `alpha_r Tr`; both
names are accepted, the records are independent).

Companions (same index): `_geometry` (geometry entity or node list),
`_element`, `_element_group`, `_element_side`, `_node`,
`_element_node` — see `condif_convection_edge_normal` for the exact
semantics (identical).

`condif_temperature` must be an unknown; 2D/3D only. The condition is
nonlinear: use `control_timestep_iterations` with enough iterations for
the Newton loop to converge (e.g. 10).

## Usage

```
condif_radiation_edge_normal <index> <alpha_r> <Tr>
condif_radiation_edge_normal_geometry <index> -<geometry_entity> <index>
```

## Example

```
condif_radiation_edge_normal_geometry 0  -geometry_line 1
condif_radiation_edge_normal 0  1.0 1.0
```

Slab with k=1, T=0 at the far end (L=1), alpha_r=1, Tr=1: steady state
T + T^4 = 1 -> T(edge) = 0.7245 (test `condif_rad`).
