# groundflow_flux_edge_normal (familia)

## Description

Distributed prescribed water flux normal to the edge of an element. The
distributed flux is translated into equivalent nodal flux on the edges of
elements (added to the pressure-dof right-hand-side).

The record `groundflow_flux_edge_normal` (index + flux value) is the main
record; `groundflow_flux_edge_normal_geometry` selects the area where it is
applied (e.g. `-geometry_line 1` in 2D), and
`groundflow_flux_edge_normal_time` optionally scales it in time.

**Attention**: this option is only available for linear and quadratic
isoparametric elements. If used INSIDE a FE mesh, the elements on each side get
the distributed flux, so the total water flux normally becomes zero (the
normals of the elements at the sides of the geometry are opposite).

## Usage

```
groundflow_flux_edge_normal <index> <flux>
groundflow_flux_edge_normal_geometry <index> <geometry_entity_name> <geometry_entity_index>
groundflow_flux_edge_normal_time <index> <time load time load ...>
```

## Main records

| Record | Meaning                                                     |
|--------|-------------------------------------------------------------|
| `groundflow_flux_edge_normal` | Distributed flux value (index + value).            |
| `groundflow_flux_edge_normal_geometry` | Area where the flux applies (2D: `-geometry_line`, 3D: `-geometry_surface`). |
| `groundflow_flux_edge_normal_time` | Time diagram scaling the flux (linear interpolation, factor 0 outside). |
| `groundflow_flux_edge_normal_sine` | Sinusoidal time load (start_time end_time freq_0 amp_0 ...). |
| `groundflow_flux_edge_normal_factor` | Polynomial in space scaling the flux (`a0 a1 ...`, a0+a1x in 1D, etc.). |

## Restriction records

| Record | Meaning                                                     |
|--------|-------------------------------------------------------------|
| `groundflow_flux_edge_normal_element` | Restrict to the listed elements.                  |
| `groundflow_flux_edge_normal_element_group` | Restrict to the listed element groups.      |
| `groundflow_flux_edge_normal_element_node` | Restrict to (element, local nodes).           |
| `groundflow_flux_edge_normal_element_node_factor` | Per-node multiplication factors for `_element_node`. |
| `groundflow_flux_edge_normal_element_side` | Restrict to (element, side) pairs.             |
| `groundflow_flux_edge_normal_node` | Restrict to the listed (global) nodes.               |

## Example

```
geometry_line 1  0. 0. 1. 0. 0.
groundflow_flux_edge_normal_geometry 0  -geometry_line 1
groundflow_flux_edge_normal 0  0.1
```

A distributed water flux of 0.1 is injected normal to the edge on line 1. In
steady state the injected flux is drained by the prescribed pressures, giving
`pres = flux * height / k` at the bottom.
