# Professional 2026 names of existing GNU records — developer

- Data part: `db_number()` aliases (`database.cc`, next to
  `global_element_dof_apply`).
- Initialization part: `input()` (`input.cc`) replaces
  `materi_history_variable` / `materi_strain_energy` by the GNU names right
  after the word is read, before `initialization_names` stores it, so the
  `check_unknown("materi_history_variables")`-style checks keep working.
- F5 T3 B: `group_materi_hyper_reduced_polynomial`,
  `group_materi_hyper_volumetric_simo_taylor`, `group_materi_umat`
  (→ `GROUP_USER_UMAT`) and `group_materi_umat_parameters`
  (→ `GROUP_USER_DATA`, the `props` of `stress_umat`) are `db_number()`
  aliases next to `group_materi_plasti_hypo_wolffersdorff`. A separate
  umat-parameters record would need its own enum and a read in
  `stress_umat`; not done because the GNU has a single user-data record.
- F5 T4 U1: `control_mesh_delete_geometry_element_group`,
  `geometry_triangle_eps_iso` and `mesh_interface_triangle_coordinates`
  are `db_number()` aliases next to
  `control_mesh_delete_geometry_move_node`. Test `alias_names_2026_3`.
- F5 T5 U1: `group_contact_spring_plasti_cohesion` and
  `group_contact_spring_plasti_friction_automatic` are `db_number()` aliases
  next to `group_contact_spring_cohesion` (`conspr.cc` already applies the
  6.605 rule: (2/3) tan(phi) of the neighbour group, 0.2 otherwise);
  `groundflow_phreatic_level_n` returns `GROUNDFLOW_PHREATICLEVEL_N` next to
  `groundflow_phreatic_level`. Test `alias_names_2026_4` (the Professional
  02-08-2026 stores the same values under the manual names).
