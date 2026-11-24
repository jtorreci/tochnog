# Professional 2026 names of existing GNU records

The manual of the Professional 02-08-2026 renames some records. The GNU
accepts the new names and stores the records under its own names (the
`.dbs` shows the GNU name):

| 2026 name (manual) | GNU record |
|---|---|
| `materi_history_variable n` (4.14, initialization part) | `materi_history_variables n` |
| `materi_strain_energy` (4.28, initialization part) | `materi_strainenergy` |
| `control_print_dof_line_node_type` (6.279) | `control_print_dof_line_method` |
| `control_print_gmsh_node_type` (6.319) | `control_print_gmsh_node_method` |
| `control_print_vtk_node_type` (6.345; -node, -node_start_refined, -node_deformed_mesh) | `control_print_vtk_node_method` |
| `element_contact_spring_direction` (6.418) | `element_contactspring_direction` |
| `element_interface_intpnt_strain_average` (6.434) | `element_interface_strain_average` |
| `element_interface_intpnt_stress_average` (6.436) | `element_interface_stress_average` |
| `post_node_rhside_ratio_dof_type` (6.951) | `post_node_rhside_ratio_unknowntypes` |
| `group_materi_hyper_reduced_polynomial` (6.678) | `group_materi_hyper_reducedpolynomial` |
| `group_materi_hyper_volumetric_simo_taylor` (6.683) | `group_materi_hyper_volumetric_simotaylor` |
| `group_materi_umat index switch` (6.757) | `group_user_umat` |
| `group_materi_umat_parameters index p0 p1 ...` (6.758) | `group_user_data` |
| `control_mesh_delete_geometry_element_group` (6.169) | `control_mesh_delete_geometry_elementgroup` |
| `geometry_triangle_eps_iso` (6.551) | `geometry_triangle_epsiso` |
| `mesh_interface_triangle_coordinates` (6.854) | `mesh_interface_triangle_coordinate` |
| `group_contact_spring_plasti_cohesion` (6.603) | `group_contactspring_cohesion` |
| `group_contact_spring_plasti_friction_automatic` (6.605) | `group_contactspring_friction_automatic` |
| `groundflow_phreatic_level_n nx ny` (6.572, no index) | `groundflow_phreaticlevel_n` |
| `groundflow_phreatic_bounda -yes/-no` (6.570, no index) | `groundflow_phreaticlevel_bounda` (see [groundflow_phreatic_level](groundflow_phreatic_level.md)) |

Tests: `validation-suite/test-2026/alias_names_2026_1`, `alias_names_2026_2`,
`alias_names_2026_3` and `alias_names_2026_4`
(build_safe checks the stored records).

## Differences with the Professional version

The `.dbs` of the GNU uses the GNU names.

`group_materi_umat_parameters` and the GNU `group_user_data` are one
record: the GNU umat takes its `props` from `group_user_data`, which the
GNU user routines `user_sigma` / `user_plasti` also read. An input that
gives both for the same group keeps the last one.

The Professional 02-08-2026 rejects the manual spelling
`mesh_interface_triangle_coordinates` (`I do not know`); it only accepts the
singular `mesh_interface_triangle_coordinate`, which the GNU also accepts.
