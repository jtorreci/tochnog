# control_print_beam_force_moment

## Implementación

- New file `print_beam_force_moment.cc`:
  `print_beam_force_moment( icontrol, task )`. Invoked from the control
  loop in `top.cc`, INSIDE the `control_print_frequency_allowed` gate
  (same as every other `control_print_*`):
  ```
  if ( frequency_allowed && db_active_index( CONTROL_PRINT_BEAM_FORCE_MOMENT, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_BEAM_FORCE_MOMENT, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
    print_beam_force_moment( icontrol, ival[0] );
  }
  ```
  `ival[0]` is the switch (`-SEPARATE_INDEX` / `-SEPARATE_SEQUENTIAL`);
  anything else is `db_error`.
- Registered in `database.cc`: `CONTROL_PRINT_BEAM_FORCE_MOMENT`
  (INTEGER, length 1), `CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES`
  (DOUBLE_PRECISION, `DATA_ITEM_SIZE` + `fixed_length = 0` so it can
  hold 4 values in 2D and 6 in 3D; `data_required` = the main record),
  `CONTROL_PRINT_BEAM_FORCE_MOMENT_SWITCH` (INTEGER, length 1,
  `data_required` = the main record). Enums (3) in `tochnog.h` /
  `tochnog-mod.h` in sync. `print_beam_force_moment.o` added to the
  makefile.

## Datos de partida (evidencia del layout)

- `ELEMENT_BEAM_MOMENT` (beam.cc:434/524): `NNOL*NDOF = 6` doubles,
  two nodes x {2 in-plane forces + 1 out-of-plane moment}, in the axes
  of the beam PLANE (the 2D mapping of `beam_3d` via `index_plane`),
  not the 3D global axes. It stores the element internal nodal force
  vector (K*u). The out-of-plane force and the in-plane moments are
  IDENTICALLY zero: the beam is a 2D element. The local frame is the
  one of `beam_2d` (local x = node1->node2 direction; local y =
  rotation of local x by +90 deg in the plane).
- `ELEMENT_TRUSS_DIRECTION` is stored with length `NDIM = 2` (in-plane)
  but the print recomputes the direction from the current NODE
  coordinates instead (robust, no stale data).
- Truss identification: element record `el[0]` is `-BEAM`, `-TRUSS` or
  `-TRUSSBEAM` (elem.cc dispatch); `ELEMENT_TRUSS_FORCE` is a scalar
  per element (positive = tension, truss.cc).

## Diseño / decisiones

- **Selection**: an element is printed when the minimum distance
  between its axis segment and the cut segment (3D segment-segment
  closest point) is < `BEAM_FORCE_MOMENT_CUT_TOL * max(1, cut_length)`
  with `BEAM_FORCE_MOMENT_CUT_TOL = 1.e-6` (hardcoded). Segments use
  the CURRENT node coordinates (deformed mesh when the analysis moves
  the mesh; with `options_mesh -fixed_in_space` the coordinates stay at
  the initial position).
- **Ordering**: the lines are sorted by ascending distance from the cut
  start (projection of the closest point on the cut onto the cut
  direction = `s * cut_length`); stable insertion sort on an index
  array (pattern of print_interface_stress). Equal distances keep
  element order.
- **Columns** (exactly the manual order): dist + 12 components in the
  local beam axes. The plane force is rotated to local axes with the
  `beam_2d` rotation matrix: `local_fx = a*f0 + b*f1`,
  `local_fy = -b*f0 + a*f1` where `(a,b)` are the components of the
  beam direction along the plane axes (`index_plane[0], index_plane[1]`
  from `group_beam_plane`, default `-x -y`; 2D fixed `(0,1)`).
  Out-of-plane components forced to zero.
- **Truss axial force**: for `-truss` and `-trussbeam` the axial
  columns are `+N` (first node) / `-N` (second node) from
  `ELEMENT_TRUSS_FORCE` — the element nodal force vector, same
  anti-symmetric pattern as the beam transverse components (decision,
  documented in the user manual). A pure `-truss` has zero beam
  moments. Elements without the required records
  (`ELEMENT_BEAM_MOMENT` for beams, `ELEMENT_TRUSS_FORCE` for trusses)
  are skipped (not yet calculated).
- **Switch**: `control_print_beam_force_moment_switch -yes` multiplies
  the 12 components by -1 (manual 6.264); `-no` or absent = no
  inversion.
- **File**: `beam_force_moment.<icontrol>` (-separate_index) or
  `beam_force_moment.<seq>` (-separate_sequential, static counter,
  pattern of print_interface_stress), opened in append mode. If no
  element crosses the cut, NO file is written (decision; a two-pass
  loop counts first).
- **Snap**: components with |v| < `1.e-6 * (1 + max|component|)` are
  snapped to zero in the output: the equilibrium solve converges to
  ~1e-6 relative, so analytic zeros (tip moment, out-of-plane
  components) would otherwise print as residual noise.
- The `_coordinates` record is mandatory: `db_error` when absent or
  with a length different from `2*ndim` or a zero-length cut.

## Detalles / gotchas

- **Segment-segment distance GOTCHA**: in the reference algorithm
  (Eberly / Real-Time Collision Detection) the first parameter `s`
  belongs to the FIRST segment and `u` to the SECOND; swapping them
  silently returns the wrong closest points (found during calibration:
  the true intersection parameter was `t=0.50083` on the cut but the
  buggy code returned `t=0.25`). The helper returns the parameter on
  the CUT segment and documents it.
- **Deformed coordinates**: the tests use `options_mesh
  -fixed_in_space` (with `options_convection -no`, required by the
  beam code) so the NODE coordinates stay at the initial position and
  the analytic distances (0.5, sqrt(1.25)) are EXACT. Without it, the
  crossing point shifts with the deformation.
- **P-delta in bmom_truss**: with an axial tension the root moment is
  F*L - N*w(L) = 0.00993, NOT F*L (second-order effect); the target
  encodes it.

## Pendiente

- 3D cut verification: the code path for `group_beam_plane` planes
  other than the default (x-y) is implemented but not covered by a
  dedicated test (the lot tests are 2D).
- The axial columns of a pure `-truss` element are `+N`/`-N` (decision);
  the Professional manual does not show a pure-truss example.
