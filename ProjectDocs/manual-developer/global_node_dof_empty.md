# global_node_dof_empty (developer)

- Enum `GLOBAL_NODE_DOF_EMPTY` (`tochnog.h` / `tochnog-mod.h`, before
  `LAST_DUMMY`), registered in `database.cc` next to `OPTIONS_ELEMENT_DOF`:
  INTEGER, length 1, `no_index`.
- `bounda.cc`:
  - `bounda_node_dormant( max_node, time_total, empty_type=0 )` gains the
    `empty_type` mode: an element is inactive when its group is not alive
    (`group_time_alive()`) or its `group_type` has only `-EMPTY` / `-NONE`
    (new static `bounda_group_type_empty()`); the
    `element_group_applied()` test is kept only for the bounda mode
    (`empty_type=0`, unchanged behaviour), and the early `NULL` return
    without `group_time` / apply records is skipped in the `empty_type` mode.
  - `node_dof_empty_apply()` (serial): reads the switch (default `-YES`),
    marks the nodes with `bounda_node_dormant(..., time_current+dtime, 1)`
    and sets every `NODE_DOF` VERSION_NEW slot to zero except the
    displacement slots (`dis_indx .. +ndim*nder`), which get the VERSION_NORMAL
    (start of step) value.
- `top.cc`: `node_dof_empty_apply()` then `dof_limit_apply()` at the end of
  the equilibrium block of every time step (after the post-solve material
  pass and the hoisted stress sync, before the step acceptance).

Verification: Professional 02-08-2026 probes f5/v5b `c_def`, `c_yes`,
`c_no` (group_type -empty), `c_tdef`, `c_tno` (dead group_time), `c2_def`,
`c2_no` (19 dofs incl. temperature, kappa, work). Tests
`global_node_dof_empty1/2/3`.

Hardcoded: the displacement exception (owner decision C101, 2026-11-03,
= Professional; the manual zeroes every dof: DIVERGENCES
GLOBAL-NODE-DOF-EMPTY-DISPLACEMENT, class A).
Pending: `ELEMENT_EMPTY` (materi_density / materi_diffusion fill) is not
counted; `group_time_fill` (birth_empty) does not exist in the GNU yet (F5
T5).
