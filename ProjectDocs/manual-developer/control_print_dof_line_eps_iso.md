# control_print_dof_line_eps_iso

## Implementación

- Stored as a CONTROL DOUBLE record (`data_length 1`). Read in
  `print_dof_line_point()` of `print_dl.cc`; default
  `DOF_LINE_DEFAULT_EPS_ISO = 1.e-3` (the Professional default, manual
  6.276); negative values -> `db_error`.
- The value is passed as the NEW optional `eps_iso` parameter of
  `point_el()` (declared `double eps_iso = 1.e-3` in `tochnog.h`).

## Diseño / decisiones

- `point_el()` uses `eps_iso` in BOTH acceptance checks of the generic
  (BAR/QUAD/HEX) branch: the isoparametric bounds
  `iso < -(1+eps) || iso > 1+eps` AND the distance check
  `dist_old > element_largest_size*eps`. The TRIA and TET branches use it
  for the barycentric bounds. `EPS_SIZE` (1.e-6, the tight
  point-to-reconstructed-position check of the TET branch and the
  element-bounding-box early reject) is NOT scaled: it is a different,
  much tighter tolerance.
- GOTCHA documented in the test (dpline_eps): with `eps_iso >= 1` a point
  clearly OUTSIDE the element can be accepted with EXTRAPOLATED shape
  functions (the distance check `dist < element_largest_size*eps` and the
  isoparametric bound both open up). This is the documented meaning of
  "increase the default value if the mesh is not exactly adjusted to the
  line" pushed to its limit.

## Detalles

- The legacy `#define EPS_ISOP 1.e-3` in point_el.cc is kept (now unused
  in that file) as documentation of the default.

## Pendiente

- None.
