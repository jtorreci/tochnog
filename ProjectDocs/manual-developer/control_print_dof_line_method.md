# control_print_dof_line_method

## Implementación

- Stored as a CONTROL INTEGER record. Read in `print_dof_line_point()`
  of `print_dl.cc`; validated against `-NODE` / `-NODE_START_REFINED`
  (default `-node_start_refined`).
- The frame is applied per node in `dof_line_interpolate()`:
  ```
  if ( method==-NODE_START_REFINED && db_active_index(
       NODE_START_REFINED, inod, VERSION_NORMAL ) )
    db( NODE_START_REFINED, inod, idum, &coords[inol*ndim], ldum, ... );
  else
    db( NODE, inod, idum, &coords[inol*ndim], ldum, ... );
  ```

## Diseño / decisiones

- The frame selects BOTH the geometry of the point-in-element test AND
  the coordinates written to the files (the line positions are in the
  same frame). This mirrors `control_print_vtk_node_method` (`-node` vs
  `-node_start_refined`; `-node_deformed_mesh` has no equivalent here).
- FALLBACK: when no `NODE_START_REFINED` record exists (geometrically
  linear analysis, no refinement/mapping) the method silently uses the
  stored `NODE` coordinates, per node. Decision: per-node fallback (a
  mixed frame can only occur in a pathological partially-refined mesh;
  the record is all-or-nothing in practice).
- In an analysis where the mesh does NOT move (`materi_displacement` +
  `-fixed_in_space` forces `options_mesh` to FIXED_IN_SPACE in top.cc)
  both frames coincide; the discriminating test uses a follow-material
  analysis (`materi_velocity` WITHOUT `materi_displacement`), where the
  mesh moves and `node_start_refined` (input as NODE-class records,
  `version_all=1`) keeps the reference.

## Detalles

- GOTCHA (test dpline_method): with `group_materi_memory -total` the
  follow-material one-element model triggers
  `matrix_inverse( old_rot )` failure ("too large distortions") in
  materi.cc; `-updated_without_rotation` runs clean and still moves the
  mesh.

## Pendiente

- None.
