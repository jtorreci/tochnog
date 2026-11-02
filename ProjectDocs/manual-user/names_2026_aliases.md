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

Test: `validation-suite/test-2026/alias_names_2026_1` (build_safe checks
the stored records).

## Differences with the Professional version

The `.dbs` of the GNU uses the GNU names.
