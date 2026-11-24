# control_print_dof_line

## Implementación

- **Output**: `print_dof_line()` in the new file `print_dl.cc` (with
  `print_dof_point()`; both share the driver
  `print_dof_line_point( icontrol, task, is_line )`). Invoked from the
  control loop in `top.cc`, inside the
  `control_print_frequency_allowed` gate (same as every other
  `control_print_*`):
  ```
  if ( frequency_allowed && db_active_index( CONTROL_PRINT_DOF_LINE, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_DOF_LINE, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
    if ( ival[0]!=-NO ) print_dof_line( icontrol, ival[0] );
  }
  ```
- **Keyword**: `control_print_dof_line` (data_class CONTROL, data_length
  1, type INTEGER) registered in `database.cc`. `ival[0]` holds the
  switch.
- **New enums**: `CONTROL_PRINT_DOF_LINE*` (9) and `CONTROL_PRINT_DOF_POINT*`
  (3) in `tochnog.h` / `tochnog-mod.h` (kept in sync), inserted between
  `CONTROL_PRINT_DOF_ID` and `CONTROL_PRINT_FREQUENCY_TIMEINTERVAL`.
- **makefile**: `print_dl.$(OBJ)` added to the object list and its rule.

## Diseño / decisiones

- **Interpolation**: a serial loop over the elements (print-time, small
  volume; no locking). For each element: optional group filter
  (`ELEMENT_GROUP` record), node coordinates of the chosen frame,
  `point_el( point, coords, weight, name, nnol, eps_iso )`. The first
  element that accepts the point wins (same semantics as
  `parallel_post_point`, with the added group filter and eps_iso). The
  interpolated value = sum of the nodal `NODE_DOF` (VERSION_NEW first,
  else NORMAL) weighted by the shape functions.
- **point_el signature change**: `point_el()` gained an optional 6th
  parameter `double eps_iso = 1.e-3` (declared in `tochnog.h`, defined in
  `point_el.cc`). The internal `EPS_ISOP` checks were replaced by the
  parameter; `EPS_SIZE` (the tight distance check of the tet branch) is
  untouched. The default equals the legacy `EPS_ISOP` so the 5 existing
  callers (post.cc, contact.cc, force.cc, map.cc) are unaffected.
- **Files**: one file per PRIMARY unknown (`db_name(dof_label[iuknwn])`,
  `iuknwn = ipuknwn*nder`) plus one per `node_dof_calcul` item
  (`post_calcul_names_without_extension[icalcul]`, the "a"+<unknown>
  convention of `print_unknowns`, e.g. `adis`, `avel`). The manual's
  "dof.index" is the dof LABEL, NOT the literal `dof.<index>` of
  `control_print_dof` (that file layout is not reused).
- **Line points**: `n` points equidistant over the TOTAL polyline length
  (decision; the manual only says "how many points will be printed"). The
  legacy `POST_LINE_N` default (5) is reused as the default of `_n`.
- **Point outside the mesh**: omitted from ALL files of that call
  (decision documented in the user manual; consistent with "the point is
  accepted to be part of an element"). No error.
- **Switch**: `-yes` behaves like `-separate_index` (same as
  `print_dof`); `-separate_sequential` uses a static counter per call
  (shared by all dof files of the call; reset per run, same limitation as
  `print_dof`'s static `dof_seq`).

## Detalles

- Records read per icontrol: `_coordinates` (required; validated
  `ncoord % ndim == 0` and at least 2 vertices), `_n` (default 5, `n<1`
  -> `db_error`), `_eps_iso` (default 1.e-3, negative -> `db_error`),
  `_method` (default `-node_start_refined`), `_element_group` (optional
  list), `_move`, `_time` (validated `-yes`/`-no`).
- `_time`: writes `# time <time_current>` once per call per file (the
  file is opened in append mode; the comment becomes the first line of
  each call's block).
- `_move`: POST_POINT_MOVE pattern (post.cc:59-66): each found point is
  displaced by the interpolated velocity times `DTIME` (VERSION_NEW
  first) and the coordinates record is PUT back. The move only runs when
  the global `materi_velocity` is initialized (manual 6.278: "should only
  be used if materi_velocity is initialised"); otherwise it is skipped
  silently.
- `nuknwn == 0` -> the routine returns without creating files.

## Pendiente

- `control_print_dof_point` has no `_element_group`/`_eps_iso`/`_method`/
  `_move` records in the Professional manual (6.281-6.283); the point
  interpolation uses the line machinery with the default settings.
- `-separate_sequential` counters are per-process statics (shared across
  icontrols, reset per run).
