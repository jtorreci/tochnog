# edge_node_restriction — `_node` / `_element_node` selection of edge loads

## What it does

The edge-load families accept a node restriction with the same index:

| family | node restriction records |
|---|---|
| `condif_heat_edge_normal` | `condif_heat_edge_normal_node`, `condif_heat_edge_normal_element_node` |
| `condif_convection_edge_normal`, `condif_radiation_edge_normal` | `..._node`, `..._element_node` |
| `groundflow_flux_edge_normal` | `groundflow_flux_edge_normal_node`, `..._element_node` |
| `force_edge`, `force_edge_normal`, `force_edge_water`, `force_edge_projected` | `..._node`, `..._element_node` |
| `support_edge_normal` | `support_edge_normal_node`, `support_edge_normal_element_node` |

Since 2026-10-06 (audit K04, Pro test `condif7`) the restriction selects
**element sides**: a side gets the load only when **all** its nodes are in
the list (for `_element_node`: the element matches and all local nodes of
the side are listed). This is the rule of the `_geometry` records ("the
total edge of an element must be inside the geometry") and the behaviour
of Tochnog Professional 02-08-2026.

## Why it matters

Before, the nodes were checked one by one, so a corner node in the list
also loaded the adjacent side that has only that node in the list. In
`condif7` (flux 1 on the left edge, nodes 1 3) the corner nodes got the
flux of the bottom and top sides as well: temperature 20 instead of the
analytic 10.

## Syntax

```
condif_heat_edge_normal          index  q
condif_heat_edge_normal_geometry index  -geometry_name geometry_index
condif_heat_edge_normal_node     index  node_0 node_1 ...
condif_heat_edge_normal_element_node index  element local_node_0 local_node_1 ...
```

- `node_i`: global node numbers. List every node of each side you want
  loaded (2 per side for linear 2D elements, 3 for quadratic ones, the
  face nodes in 3D).
- `element`, `local_node_i`: one element and its local node numbers.
- A single node selects no side (no load).

## Example

```
condif_heat_edge_normal      0  1.0
condif_heat_edge_normal_geometry 0  -world
condif_heat_edge_normal_node 0  1 3 4
```

On a 10 x 10 quad4 (nodes 1 2 / 3 4, temp 0 at nodes 2 and 4) this loads
the left side (1-3) and the top side (3-4), not the bottom side (1-2):
temp(1) = 12, temp(3) = 18, as the Professional. Test:
`validation-suite/test-2026/edge_node_restriction_heat1.dat`
(`edge_node_restriction_force1.dat` for `force_edge_node`).

## Differences with the Professional version

None known. Measured equal for heat flux, groundflow flux, convection,
`force_edge`, `force_edge_normal` and `support_edge_normal`; radiation,
water and projected edge loads share the same code (not probed).
