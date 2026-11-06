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
