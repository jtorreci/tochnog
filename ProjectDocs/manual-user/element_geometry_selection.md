# Element selection by a geometry

`geometry_method`, `control_mesh_delete_geometry_direct`,
`control_mesh_delete_geometry_method`,
`control_mesh_delete_geometry_projection_type`,
`control_mesh_refine_locally_minimal_size`,
`control_mesh_refine_locally_not_method`,
`control_mesh_refine_locally_only_method` (manual Professional 6.532,
6.167, 6.171, 6.173, 6.224, 6.226, 6.228).

Several commands select elements with a geometry entity. These records
choose when an element counts as "in" the geometry, and adjust the delete
and local refinement commands.

## Selection methods

| method | the element is selected when |
|---|---|
| `-all` | all its nodes are in the geometry (completely inside) |
| `-any` | at least one of its nodes is in the geometry (partially inside) |
| `-average` | its middle coordinate (the average of its node coordinates) is in the geometry |

## Records

| record | meaning | default |
|---|---|---|
| `geometry_method index method` | method `-all`, `-any` or `-average` for the geometry with the same index, used by `element_geometry_present` (`print_element_geometry_present -yes`) and `dependency_geometry` | `-all` |
| `control_mesh_delete_geometry_method index method` | method `-all`, `-any` or `-average` of `control_mesh_delete_geometry` with the same index | `-all` |
| `control_mesh_delete_geometry_projection_type index type` | `-project_inside`: everything inside the geometry is deleted (for a circle: the whole disk plus the tolerance); `-project_exact`: only what lies within the tolerance of the geometry edge | `-project_inside` (see the differences below) |
| `control_mesh_delete_geometry_direct index switch` | `-yes`: the elements are deleted at once, in the first step, instead of slowly over the timesteps; a `control_mesh_delete_geometry_factor` is then not used | `-no` |
| `control_mesh_refine_locally_not_method index method` | `-all` or `-any`: an element is excluded by `control_mesh_refine_locally_not` when all / any of its nodes are in that geometry | `-all` |
| `control_mesh_refine_locally_only_method index method` | `-all` or `-any`: an element may be refined by `control_mesh_refine_locally_only` when all / any of its nodes are in that geometry | `-all` |
| `control_mesh_refine_locally_minimal_size index minimal_size` | an element whose size (the largest distance between two of its nodes) is below `minimal_size` is not refined | 0 |

`geometry_method` does not change `control_mesh_delete_geometry`,
`control_mesh_refine_locally_not` / `_only` or `area_element_group`: these
commands have their own method records (as in the Professional).

A per-geometry `geometry_projection_type` (6.543) still overrides
`control_mesh_delete_geometry_projection_type` for that geometry.

With `-average` the middle coordinate is tested as a point: the node-only
geometry options (`geometry_boundary`, `geometry_mpc`,
`geometry_element_geometry`, `geometry_element_group`) and
`geometry_node_type` do not apply to it.

A value outside the lists above is an input error (`Error detected for
data item : <record>`).

## Example

```
geometry_quadrilateral 7 0.4 -0.01 1.01 -0.01 0.4 1.01 1.01 1.01 0.001
geometry_method 7 -any
print_element_geometry_present -yes

control_timestep 10 0.1 0.4
control_mesh_delete_geometry 10 -geometry_quadrilateral 7
control_mesh_delete_geometry_method 10 -average
control_mesh_delete_geometry_direct 10 -yes

control_mesh_refine_locally 20 100.
control_mesh_refine_locally_only 20 -geometry_quadrilateral 7
control_mesh_refine_locally_only_method 20 -any
control_mesh_refine_locally_minimal_size 20 0.1
```

Tests (validation-suite/test-2026): `geometry_method1`-`5`,
`geometry_method_error1`, `delete_geometry_method1`-`3`,
`delete_geometry_method_error1`, `delete_geometry_projection1`-`3`,
`delete_geometry_projection_error1`, `delete_geometry_direct1`-`2`,
`delete_geometry_direct_error1`, `refine_locally_only_method1`-`2`,
`refine_locally_not_method1`-`2`, `refine_locally_method_error1`,
`refine_locally_minimal_size1`-`2`.

## Differences with the Professional version

- `control_mesh_delete_geometry` with a circle: the GNU deletes the
  filled circle by default (manual 6.166: "elements in the total inner area
  of the circle radius plus its tolerance will be deleted") and follows
  `control_mesh_delete_geometry_projection_type`. The manual 6.173 gives
  `-project_exact` as the default, and the Professional 02-08-2026 always
  tests the circle edge only, ignoring the record and
  `geometry_projection_type` (DIVERGENCES DELETE-PROJECTION-TYPE; the
  default is an open owner decision, BACKLOG C123).
- Invalid method / switch values are input errors in the GNU; the
  Professional ignores them, runs a refine `-average` as `-any`, and stops
  with an assertion for `geometry_method -only` (DIVERGENCES
  U4-METHOD-VALUES).
- The GNU does not print `element_delete_factor` in the `.dbs`; the
  Professional does (the slow deletion itself is the same, BACKLOG C124).
- Before F5 T4 U4 the GNU used `-any` for `control_mesh_refine_locally_only`;
  it now uses `-all` like the manual and the Professional.
