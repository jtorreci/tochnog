# force_edge_diagram

## Description

`force_edge_diagram index vertical_displacement_0 factor_0
vertical_displacement_1 factor_1 ...` (manual 6.455): multiplication
factors for the `force_edge` record with the same index, as a function of
the vertical displacement (`disy` in 2D, `disz` in 3D). Give the
displacements from low to high; in geotechnical models they are usually
negative.

## Effect on the calculation

- The factor is evaluated per edge node, from that node's vertical
  displacement, and multiplies the node's share of the edge force (as
  `force_edge_node_factor` does).
- The displacement is the one of the previous equilibrium iteration, so the
  load follows the settlement within the iterations of a step and from step
  to step (Professional 02-08-2026).
- Outside the diagram the factor is 0, like `force_edge_time` outside its
  times: cover every displacement the model can reach.
- It multiplies the other factors (`force_edge_factor`, `force_edge_time`,
  `force_edge_node_factor`, ...).
- It needs `materi_displacement`; without it GNU stops with `Error:
  force_edge_diagram needs materi_displacement`.

## Example

```
force_edge 0  0. -1.
force_edge_geometry 0  -geometry_line 1
force_edge_diagram 0  -1. 0.  0. 1.
```

Here the load drops linearly from 1 at no settlement to 0 at a settlement
of 1.

Self-checking inputs in `validation-suite/test-2026`: `force_edge_diagram1`,
`force_edge_diagram2` (outside the diagram), `force_edge_diagram_error1`.
