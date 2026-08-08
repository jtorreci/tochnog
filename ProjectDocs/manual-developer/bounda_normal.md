# bounda_normal

## Files and functions

- `bounda.cc` — within `bounda()`.
  - Local `bounda_normal_vec[3]` declared at line 48.
  - Normal vector read per iboun: `db( BOUNDA_NORMAL, iboun, idum,
    bounda_normal_vec, ldum, VERSION_NORMAL, GET_IF_EXISTS )` (lines 200-202).
  - Projection applied at the end of the node dof loop (lines 560-579), after
    all velocity dofs of the node have been fixed.
- `database.cc` — keyword registration (lines 212-216): name
  `bounda_normal`, `type = DOUBLE_PRECISION`, `data_length = 3`,
  `data_class = BOUNDA`, `data_required = BOUNDA_UNKNOWN`.
- `tochnog.h` — enum `BOUNDA_NORMAL` (line 145).
- `tochnog-mod.h` — mirror enum `BOUNDA_NORMAL` (line 138), must stay in sync.

## Implementation details

- The normal vector is read with `GET_IF_EXISTS`, so the projection is
  inactive when no `bounda_normal` record exists for `iboun`.
- The projection is applied only when `materi_velocity` is set and the vector
  is non-zero (bounda.cc:562-563).
- The velocity vector `v` (components `new_node_dof[vel_indx+idim*nder]` for
  `idim=0..ndim-1`) is projected onto the plane by removing the normal
  component: `vn = v·n` then `v -= (vn/|n|²)·n`. This is a per-node operation.
- The projection runs at the END of the dof loop of the node (after the
  prescribed dofs and derived `MATERI_DISPLACEMENT`/
  `MATERI_VELOCITY_INTEGRATED` velocities are written), so the full velocity
  vector is available before the plane projection.
- Uses `ndim`, so in 2D only two components are projected and in 1D only one.

## External dependencies

None. Internal database API only; reuses `materi_velocity`, the velocity dof
offset `vel_indx` and the global `ndim`.

## Hardcoded parameters / pending refactorings

- The normal vector length is hardcoded to 3 components; only the first `ndim`
  are used.
- The zero-vector check and the `nn>0.` guard handle a zero/undefined normal
  (the node then stays unconstrained), but the vector is not normalized once
  and reused; the normalization factor `nn` is recomputed per node.
- Only the velocity dofs are projected; other physics that carries the velocity
  (e.g. Maxwell) is not adjusted here.
