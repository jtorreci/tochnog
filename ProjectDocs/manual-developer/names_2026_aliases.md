# Professional 2026 names of existing GNU records — developer

- Data part: `db_number()` aliases (`database.cc`, next to
  `global_element_dof_apply`).
- Initialization part: `input()` (`input.cc`) replaces
  `materi_history_variable` / `materi_strain_energy` by the GNU names right
  after the word is read, before `initialization_names` stores it, so the
  `check_unknown("materi_history_variables")`-style checks keep working.
