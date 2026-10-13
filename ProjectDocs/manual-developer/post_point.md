# post_point (developer)

## Where

- `post.cc`, `post()`: loops over the active `POST_POINT` records, calls
  `parallel_post_point()` through `parallel_sys_routine()` and stores
  `POST_POINT_DOF`.
- `post.cc`, `parallel_post_point()`: for each active element, reads the
  node coordinates from `NODE_START_REFINED` and calls `point_el()`; the
  first element that contains the point wins; `post_point_dof` = sum of
  `weight[inol] * node_dof` (VERSION_NEW when present).
- `point_el.cc`, `point_el()`: point location and shape-function weights
  per element type:
  - tria3/tria6: barycentric projection (closed form);
  - tet4/tet10: volume coordinates (closed form);
  - prism15: finite-difference descent over (L1, L2, zeta);
  - bar/quad/hex (2, 3, 4 nodes per direction): finite-difference
    descent over the isoparametric coordinates with the
    `interpolation_polynomial()` Lagrange bases (xi fastest, zeta
    slowest), full isoparametric geometry.

## Implementation details

- The descent stops when every iso step is below `EPS_SIZE` (1e-6); its
  final accepted step leaves the point about 1e-8 iso units short
  (`mpc.cc` polishes its tie weights with a Gauss-Newton step for this
  reason; `post_point` does not).
- Acceptance: iso coordinates within `[-1-eps_iso, 1+eps_iso]` (default
  `eps_iso` 1e-3, `tochnog.h`) and a residual distance below
  `eps_iso * element size`.
- `point_el()` reads coordinates from the caller; every caller passes
  `NODE_START_REFINED`, so a mesh generator that leaves that record stale
  breaks the location (F4 T6f C57: the quad9 -> hex27 extrusion, fixed in
  `mesh.cc` `mesh_extrude()`).

## Hardcoded parameters

- `EPS_ISOP` 1e-3, `EPS_SIZE` 1e-6, `MITER` 10000, `DELTA` 1e-4
  (`point_el.cc`).

## Differences with the Professional version

- C78, deliberate divergence (project owner 2026-10-13): the Professional locates points in quad9 / hex27 with the
  corner nodes only (sub-parametric location, measured 2026-10-13 on a
  curved quad9 and a distorted hex27: iso within about 3e-6 of the exact
  bilinear / trilinear inversion), then interpolates with the quadratic
  shape functions. The GNU keeps the exact isoparametric location
  (`point_el()` unchanged); a post point on a node returns the nodal
  value.
- `tet10`: quadratic interpolation (Professional linear, TET10-PPOINT).

## Pending refactorings

- An analytic-Jacobian Newton polish (as in `mpc.cc`) would remove the
  1e-8 shortfall of the descent (prism15 `post_point` 1.000000261 vs
  1.0).
