# control_print_dof_point

## Implementación

- **Output**: `print_dof_point()` in `print_dl.cc`, which calls the
  shared driver `print_dof_line_point( icontrol, task, is_line=0 )`
  with `npoints = 1` and the point coordinates record. Invoked from the
  control loop in `top.cc` inside the `control_print_frequency_allowed`
  gate (same as `control_print_dof_line`).
- **Keyword**: `control_print_dof_point` (data_class CONTROL,
  data_length 1, type INTEGER) registered in `database.cc`.

## Diseño / decisiones

- The point family shares ALL the line machinery (interpolation with
  group filter/eps_iso/method, file naming, `# time` comment). The
  Professional manual defines no `_element_group`/`_eps_iso`/`_method`/
  `_move` records for the point (6.281-6.283), so the defaults apply:
  all groups, eps_iso 1.e-3, `-node_start_refined`, no move.
- A point outside the mesh: no data line for that call (same omission
  decision as the line family).
- `node_dof_calcul`: interpolated like `node_dof` (same weights) and
  written to files named `post_calcul_names_without_extension[icalcul]`
  (the print_unknowns convention, e.g. `adis`, `avel`). This mirrors the
  manual statement that dof_point prints "values of the node_dof records
  and node_dof_calcul records" (6.281) — note that the legacy
  `print_dof()` of `control_print_dof` does NOT print `node_dof_calcul`
  (documented difference; the line/point family does, following the
  manual).

## Detalles

- The point coordinates record is validated (`ncoord >= ndim`).
- One file per PRIMARY unknown (`db_name(dof_label[iuknwn])`), so a
  2D elastic model writes `velx.<ext>`, `vely.<ext>`, `disx.<ext>`,
  `disy.<ext>` (the manual example `temp.10`/`velx.10`).
- GOTCHA of the calcul test (dpoint1): with `materi_displacement` +
  `-fixed_in_space` the mesh does not move (F = I), so a
  `-materi_strain_total -average` calcul is exactly 0; the test uses
  `-materi_displacement -average` instead (mean of the displacement
  vector = (disx+disy)/2, non-zero and time-varying).

## Pendiente

- None.
