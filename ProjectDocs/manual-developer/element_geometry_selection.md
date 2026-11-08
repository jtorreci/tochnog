# Element selection by a geometry (developer)

Records `geometry_method`, `control_mesh_delete_geometry_direct`,
`_method`, `_projection_type`, `control_mesh_refine_locally_minimal_size`,
`_not_method`, `_only_method` (manual Professional 6.532, 6.167, 6.171,
6.173, 6.224, 6.226, 6.228; F5 T4 U4).

## Where

- `tochnog.h` / `tochnog-mod.h`: seven enums right after
  `GEOMETRY_TETRAHEDRAL` (same order in both headers); prototypes of
  `element_in_geometry()` and `geometry_method_get()` after
  `element_geometry_present_calculate()`.
- `database.cc`: `GEOMETRY_METHOD` (class GEOMETRY, INTEGER 1) after
  `GEOMETRY_TETRAHEDRAL`; the delete records after
  `CONTROL_MESH_DELETE_GEOMETRY_MOVENODES` (required
  `CONTROL_MESH_DELETE_GEOMETRY`); the refine records after
  `CONTROL_MESH_REFINE_LOCALLY_ONLY` (required `CONTROL_MESH_REFINE_LOCALLY`,
  `_NOT`, `_ONLY`); the keyword value `project_exact` (enum `PROJECT_EXACT`
  existed, unnamed) next to `project_inside`.
- `geometry.cc`:
  - `geometry_method_get(geometry_entity)`: reads `GEOMETRY_METHOD` at the
    geometry index (default `-ALL`), `db_error` for another value;
  - `element_geometry_coord()` (static): average node coordinate with the
    `node_type` rules of `geometry()` (`PLUS_DISPLACEMENT` adds the
    displacements);
  - `element_in_geometry(element, geometry_entity, method, node_type,
    projection_type, node_in_geometry, version)`: the shared selection;
    `-ALL` / `-ANY` per node (from the caller's `node_in_geometry` array
    when given, else `geometry()` per node), `-AVERAGE` one coordinate-only
    `geometry()` call (inod -1);
  - `element_geometry_present_calculate()` uses it with
    `geometry_method_get()` of each geometry.
- `group.cc`, `dependency_geometry_contains()`: same helper with the
  `geometry_method` of the dependency geometry.
- `delete.cc`, `delete_geom()`: reads `_direct` (`-yes` skips the
  time-interpolated `element_delete_factor`, so the element is removed in
  the first step), `_method` (passed to `element_in_geometry` with the
  precomputed `node_in_geometry` array) and `_projection_type` (maps to
  `PROJECT_INSIDE` / `PROJECT_EXACT`; without the record the projection
  type stays `CONTROL_MESH_DELETE_GEOMETRY`, the filled test of
  `geometry()`). The stop records keep their all-nodes test.
- `refine_l.cc`, `refine_locally()`: `_not_method` / `_only_method` through
  `element_in_geometry()` (`PROJECT_EXACT`, `NODE_START_REFINED`), replacing
  the hard-coded not = all nodes / only = any node; `_minimal_size` marks an
  element as excluded (`all_not`) when its largest node distance (`NODE`
  coordinates) is below the value. Only the refining pass uses them; the
  correcting pass is unchanged.

## Measured on the Professional 02-08-2026

Probes in a 4x4 quad4 mesh (quadrilateral x 0.4..1.01 over element 1, circle
centre 2 2 radius 1.45 / 0.5), two quad4 with a dependency point, eight
tria3 refined with `-nothing`:

- `element_geometry_present`: `-any` 1 2 5 6, `-average` 1, default none;
- `dependency_geometry`: `-any` 300 / 300, `-average` 300 / 2000, default
  1000 / 2000;
- `geometry_method` is ignored by `control_mesh_delete_geometry`,
  `control_mesh_refine_locally_not` / `_only` and `area_element_group`;
- delete `-any` keeps 3 4 7 8 9-16, `-average` deletes element 1;
- delete with a circle: nothing deleted with any projection record or
  `geometry_projection_type` (edge test only); `-any` deletes everything
  (every element has a node on the edge);
- `_direct -yes`: elements removed in the first step, also with
  `_factor 0.2 0.6`; `-no`: element_delete_factor 0.75 / 0.5 / 0.25 then
  removal (factor 0.2 0.6: 0.7 / 0.6 / 0.5 / 0.4, kept);
- refine: default only = `-all` (16 elements / 14 nodes), only `-any` 32 / 25,
  not default and `-all` 28 / 22, not `-any` 8 / 9, minimal size 1.5 8 / 9,
  1.4142 and 1.0 32 / 25;
- invalid values: refine `-average` behaves as `-any`, delete method /
  projection / direct values are ignored, `geometry_method -only` stops with
  "Program error, assertion failed".

## Hardcoded

- The middle coordinate is the arithmetic mean of the element nodes (equal
  to the iso-parametric centre for linear elements, not for quadratic
  serendipity elements).
- Element size for `_minimal_size`: largest distance between any two nodes
  (`NODE` coordinates), compared with `<`.

## Pending

- Default of `control_mesh_delete_geometry_projection_type` (6.166 vs
  6.173 / Professional): owner decision C123 (DIVERGENCES
  DELETE-PROJECTION-TYPE).
- `element_delete_factor` is not printed by the GNU (`external = 0`), the
  Professional prints it: C124.
- `-average` ignores the node-only filters and `geometry_node_type` of the
  geometry (the middle coordinate is not a node).
- Other all-nodes selections (`control_mesh_split_only`,
  `control_data_initeldof_geometry`, `control_mesh_cut_geometry`) keep their
  own test; they could use `element_in_geometry()`.
