# global_element_dof_apply — developer

- `db_number()` alias (`database.cc`) of `OPTIONS_ELEMENT_DOF` (INTEGER,
  `no_index`), global `options_element_dof` (initia.cc, default `-YES`).
- New switch word `-only` (enum `ONLY`, after `LAMBDA`). `top()`
  (`top.cc`) turns `-ONLY` into `element_dof_only = 1` (initia.cc) and
  `options_element_dof = -YES`, so every `options_element_dof==-YES`
  consumer keeps storing `ELEMENT_INTPNT_DOF`.
- `element_dof_only` skips the node update of non-principal dofs in
  `parallel_new_dof_diagonal()` (`dof.cc`; the displacement is still
  integrated from the velocity afterwards) and the hoisted post-solve
  nodal stress assignment in `top()`.
- Not changed: the derivative averaging of `element_loop()` (`elem.cc`).

## Pending

- BACKLOG C93: `post_element_force` with `-no` (no `ELEMENT_INTPNT_DOF`) gives 0.
