# force_edge_projected

## Description

Terzaghi projected tunnel load (manual Professional 6.478): a linear
ground stress field is projected on the excavation boundary and applied
as the release load. The master record (Professional syntax; the parser
maps it to `force_element_edge_projected` internally):

```
force_edge_projected <index>
    ph(0,0,0) ph_grad_x ph_grad_y [ph_grad_z]
    pv(0,0,0) pv_grad_x pv_grad_y [pv_grad_z]
    factor_normal factor_tangential
    vertical_dir_downward_x vertical_dir_downward_y [vertical_dir_downward_z]
    [tunnel_dir_x tunnel_dir_y tunnel_dir_z]        (3D only)
```

(10 values in 2D, 16 in 3D.)

Physics: `ph` is the horizontal ground stress (perpendicular to the
tunnel axis), `pv` the vertical one — both linear in space through
their gradients. With the edge outward normal n, tangent t, the
downward unit vector vd and hd = tunnel x vd (in 2D hd is the
in-plane horizontal):

- radial stress: `sig_radial = ph (n.hd)^2 + pv (n.vd)^2`
- tangential shear: `sig_tangential = ph (t.hd)(n.hd) + pv (t.vd)(n.vd)`

applied as `factor_normal*sig_radial*n + factor_tangential*sig_tangential*t`
(the load pushes the boundary outward into the excavation).
Horizontal stress along the tunnel axis is not included.

Companions (same index): `_geometry` (required; the excavation edge),
`_time`/`_sine` (temporal factor), `_factor` (spatial polynomial),
`_element`, `_element_group`, `_element_side`, `_node`,
`_element_node`, `_node_factor` — same semantics as the other
force_edge families.

Attention (manual): for linear/quadratic isoparametric elements; used
INSIDE a mesh the contributions of both sides cancel.

## Example

```
force_edge_projected_geometry 1  -geometry_line 1
force_edge_projected 1
    10. 0. 0.
    20. 0. 0.
    1.0 0.0
    0. -1.
```

ph=10, pv=20, full normal projection. On a VERTICAL wall the radial
stress is ph=10; on a HORIZONTAL wall pv=20 (test `fproj_tunnel`:
sigxx = 10 EXACT, sigyy = 20 — the 2:1 ratio proves the projection per
wall orientation; a uniform pressure would give 1:1).
