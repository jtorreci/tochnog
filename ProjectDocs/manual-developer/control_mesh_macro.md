# control_mesh_macro (tensor macros numbering) — developer notes (F4 T8)

## Where it is implemented

`macro.cc` `macro()`, branch `macro_item==-BAR || -RECTANGLE || -BRICK`:
the node loop and the element loop.

## Implementation details

- Per-dimension counts `ncount[]` and lengths `lcount[]` (1D: the `nz` /
  `lenz` of the old code, 2D: `ny,nz`, 3D: `nx,ny,nz`).
- Nodes: loops `iz` (slowest), `iy`, `ix` (fastest); node offset
  `ix + n_x*(iy + n_y*iz)` from `inod_first`.
- Elements: same loop order; local nodes in tensor order (x fastest):
  quad4 `n0, n0+1, n0+n_x, n0+n_x+1`, hex8 the same plus `+n_x*n_y`.
- The old loops ran x slowest (the transposed numbering, and local node
  1 of a quad4 along y). Professional 02-08-2026 probes (F4 T8: 1D bar,
  2D 3x3 rectangle, 3D 3x2x2 brick) give x fastest; the GNU nodes and
  elements are now identical to the Professional `.dbs`.

## External features / libraries

None.

## Hardcoded parameters

None.

## Pending refactorings

- `node_remesh_allowed` keeps the old loop-role mapping (3D: boundary in
  direction d -> slot d; 2D: x boundary -> slot 1 and the y boundary is
  dropped; 1D: never). It looks like a legacy bug (2D slot 0 is never
  set) and is used by `remesh.cc`; left unchanged to avoid a remesh
  behaviour change without a Professional reference.
- The other macro items (`-circle_hollow`, `-cylinder_hollow`, geometry
  macros) were not compared with the Professional numbering.
