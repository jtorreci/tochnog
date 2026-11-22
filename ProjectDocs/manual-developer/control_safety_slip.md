# control_safety_slip — developer

## Files and functions

- `safety_slip.cc` (new module, F5 T6): `safety_slip()` called from
  `step_close()` (`top.cc`, after `post()`, converged steps only, before
  control_print and GiD).
  - `safety_slip_add_piece()`: points, normals (direction turned
    counterclockwise), trapezoid weights and point pairs of a line piece.
  - `safety_slip_evaluate()`: serial loop over all elements per point
    (`point_el`, eps 1e-3, current `NODE` coordinates); element nodal
    stresses `element_dof` interpolated and averaged over the containing
    elements; c, phi via `get_group_data` (mohr_coul, mohr_coul_direct;
    phi in radians); sigma_n > 0 points excluded.
  - `safety_slip_zone_lines()`, `safety_slip_zone_circles()`: one zone
    per type and index (first/last values, numbers of variations).
  - `safety_slip_variation()`: the points of one variation of a zone.
  - `safety_slip_shrink()`: the zone of the next refinement pass.
  - `safety_slip()`: passes over all zones (types in alphabetical order),
    results per index, critical flags per index or per safety_slip_set.
  - `safety_slip_add_circle()`: circle points.
  - `safety_slip_put()`: the internal `safety_slip_surface_*` records.
- `database.cc`: `CONTROL_SAFETY_SLIP` (INTEGER, CONTROL),
  `SAFETY_SLIP_MULTI_LINEAR*`, `SAFETY_SLIP_COMBINED_LINEAR*`,
  `SAFETY_SLIP_CIRCLE_GRID_*`, `SAFETY_SLIP_CIRCLE_LINE_*`,
  `SAFETY_SLIP_SET`, `SAFETY_SLIP_SET_RESULT`,
  `SAFETY_SLIP_SURFACE_*` (length
  `SAFETY_SLIP_MPOINT` = 2000 points, `tochnog.h`).
- `print_g6.cc`: `safety_slip_print_gid_mesh()` at the end of the mesh
  block (first print or -separate), `safety_slip_print_gid_results()`
  after the nodal results (not for eigen analyses); both in
  `safety_slip.cc`. The GiD 5 writer (`print_g5.cc`) is not extended.
- `makefile`, `makefile - Windows-Mingw64`: `safety_slip` object.

## Implementation details

Measured on the Professional 02-08-2026 (odd/tasks/f5-safety-slip.md,
L-T6-measure): uniform strain fields on a 4x4 quad4 mesh, results exact;
the point stresses equal the average of the element-wise interpolated
`element_dof` to 1e-12; tension and outside points excluded; surface
numbering over all surfaces of the step. The per-point limit is the
average of the element limits (equal to c - avg(sigma_n) tan phi inside
one group). Plastic Mohr-Coulomb states differ between the two binaries
(a known MC difference), so tests use elastic states.
