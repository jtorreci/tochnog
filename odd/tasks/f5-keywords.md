# F5 — missing Professional keywords

Objective: implement (or deliberately discard) the 257 Professional keywords the GNU parser
still rejects. Branch `dev/f5-keywords` (from `dev/masin-kernel-surface` `210e6f6`); delivery
`exception-ok` (local branch, merged in parts into `dev/masin-kernel-surface` with owner
approval, as F4); runner `TN_SKIP_SUITE=1 scripts/build_safe.sh --clean` then
`TN_CHECK_ONLY=1 scripts/build_safe.sh`.

## Specs
- S1 Scope: the missing list `ProjectDocs/keyword-gap-2026-11-01.txt` minus the 63 GiD names
  (discarded): 257 keywords, re-extracted from the manual of the Professional 02-08-2026
  (`doc/user.pdf`). Plan item (owner, 2026-10-10 plan): "F5: ~258 missing keywords".
- S2 Rule order (owner): the manual and analytic results decide; the Professional binary is the
  oracle only when the manual is silent; every divergence gets an owner decision and is
  documented (DIVERGENCES-PROFESSIONAL.md, SEGUIMIENTO-CONVERGENCIA.md, the manuals'
  "Differences with the Professional version", BACKLOG-OPEN-ISSUES.md).
- S3 Verification (AGENTS.md): before implementing a keyword, search a usage example (Professional
  manual, Professional repo/Drive, sfnet suite); verify against it, otherwise build a dedicated
  test; record how it was verified.
- S4 Documentation (AGENTS.md): every keyword gets `ProjectDocs/manual-user/<name>.md` and
  `ProjectDocs/manual-developer/<name>.md` (a family may share one page per the existing
  convention), both README indexes, and its SEGUIMIENTO line/row; a discarded keyword goes to the
  "Descarte de features" table with the reason.
- S5 Every implemented keyword is accepted by the probe (`scripts/keyword_probe.cc`) and the
  KEYWORD-GAP-REPORT counts are refreshed at the end of each wave.

## Tasks
- [x] T0 Refresh the gap against the 02-08-2026 manual (S1, S5). Route: inline. 1156 names,
  824 accepted + 6 literals, 320 missing = 63 GiD + 257. List `keyword-gap-2026-11-01.txt`,
  report update in KEYWORD-GAP-REPORT.md.
- [ ] T1 W1 switches and control (35): control_check_data, control_element_group,
  control_element_group_apply, control_input, control_mesh_element_group_apply, control_post,
  control_post_apply, control_post_element_force, control_reset_value_logarithmic_first,
  control_slide_damping_apply, control_system_call, control_timestep_adjust_minimum_iterations,
  control_timestep_iterations_automatic_minimum_maximum_wished, control_timestep_until_data,
  control_timestep_until_maximum, control_timestep_until_minimum, control_truss_rope_apply,
  control_zip, element_group_apply, global_element_dof_apply, groundflow_apply,
  group_plasti_apply, interface_gap_apply, linear_calculation_apply, materi_plasti_tension_apply,
  materi_plasti_visco_apply, mesh_element_group_apply, message, mrange, mstring,
  print_data_name, print_precision, processors_maximum, processors_partition, truss_rope_apply.
  Route: inline (per-unit commits); map by a read-only explorer (mapping trigger: 35 keywords
  across ~20 files). control_slide_damping_apply moves to T5 (switch of slide_damping). Units:
  - [x] U1 apply switches: materi_plasti_tension_apply, materi_plasti_visco_apply,
    group_plasti_apply, truss_rope_apply, control_truss_rope_apply, interface_gap_apply. Done `467cdbd`.
  - [ ] U2 aliases and run-level switches: control_reset_value_logarithmic_first (+ data.cc
    formula `log(b+c*x*d+e*x)` vs manual `ln((b+cx)/(d+ex))`: owner decision),
    global_element_dof_apply (`-only` new), message, print_data_name, control_check_data,
    control_zip, control_system_call (runs a user script: owner decision), processors_maximum
    (Professional default -yes vs GNU 1 thread: owner decision).
  - [x] U3 post gates: control_post, control_post_apply (+ global post_apply, registered, unread),
    control_post_element_force. Done `ea96bf6`.
  - [ ] U4 timestep loop: control_timestep_iterations_automatic_minimum_maximum_wished (GNU
    hardcodes 2/20/4, manual default 2 8 4), control_timestep_adjust_minimum_iterations,
    control_timestep_until_data / _maximum / _minimum.
  - [ ] U5 element-group selection: control_element_group, element_group_apply,
    control_element_group_apply, mesh_element_group_apply, control_mesh_element_group_apply,
    processors_partition, print_precision.
  - [ ] U6 input / database rewrite: mstring, mrange, control_input, groundflow_apply,
    linear_calculation_apply.
- [ ] T2 W2 output records (69): `element_*`, `post_*`, `print_*`, `control_print_*_node_type`,
  `materi_strain_*` / `materi_history_variable` initial dofs, `dof_*`, `global_*`,
  node_bounded_index (list in the 2026-11-01 file). Route: inline.
- [ ] T3 W3 materials (27): `group_materi_*` and `materi_plasti_*` (umat, stress_null, hypo
  niemunis/wolffersdorff variants, borja_tamagnini, elasti_c, hyper reduced_polynomial /
  simo_taylor, visco_power, residual_stiffness, mpc, diprisco_density, history_variable_user,
  phimob, maximum_iterations, hypo_substepping). Route: per family.
- [ ] T4 W5 geometry, mesh tools, mpc, dependency (69). Route: per family.
- [ ] T5 W6 other physics (25): baseline correction, force_edge_diagram, force_gravity_geometry,
  force_volume_element_group_0, groundflow_phreatic_*, group_beam_*, group_contact_spring_*,
  group_dof_initial*, group_time_fill, group_truss_elasti_elongation_force_diagram,
  slide_damping, slide_user, strain_volume_*. Route: per family.
- [ ] T6 W7 safety_slip (32): self-contained slip-surface safety-factor feature, planned
  separately. Route: to decide.

## Log
- L1 2026-11-01 owner, after "Lo próximo es F5, los ~258 keywords que faltan.": "vengaaaa!"
- 2026-11-01 T0: probe rebuilt from the `7da9e21` objects; the 2026-09-26 list gave 321 still
  missing (only `group_integration_method_reduced_factor` added since). The 02-08-2026 manual
  renames `*_node_method` to `*_node_type` (print_vtk/gmsh/gid, control_print_*, new
  control_print_dof_line_node_type), `global_post_point_node_type` to `post_point_node_type`
  and drops `mesh_delete_geometry_moving`. Next: T1.
- 2026-11-02 T1 map (read-only explorer): 6 units U1-U6 (see T1); only control_input,
  control_post, control_timestep_iterations_automatic_minimum_maximum_wished and
  linear_calculation_apply have Professional usage examples; owner decisions pending for U2
  (logarithmic reset formula, control_system_call, processors_maximum default).
- 2026-11-02 U1 (route inline). Professional probes (29 variants): tension / visco / group
  switches = Professional, control overrules global; the control tension form now also gates
  group_materi_plasti_tension (= Professional, GNU before ignored it). The Professional keeps the
  rope with truss_rope_apply -no and rejects interface_gap_apply: GNU follows the manual
  (DIVERGENCES ROPE-APPLY, IFACE-GAP-APPLY, class D). RED: 12 new tests rc 1 on 4b80565. Differential
  vs 4b80565 (75 inputs): only new keywords and time_calculation lines. build_safe 364 OK / 0
  FALLO, 226/226; audit_db_buffers rc 0. Commit 467cdbd. Next: U3 (post gates; U2 waits for owner
  decisions).
- 2026-11-02 U3 (route inline). Professional probes with control_print_history series per
  control index (11 variants): post_apply -no creates no post record; control_post_apply -no
  freezes them in its index (control overrules global); control_post_element_force -no freezes
  only the force; an index without timesteps evaluates the post records (control_post -yes/-no
  alike); the Professional re-evaluates the post records at the end of the run. GNU: gate in
  step_close, stale flags, end-of-run re-evaluation of stale records only (runs without the
  switches byte-identical). 5 tests + 5 history/record checks in build_safe; build_safe 374 OK /
  0 FALLO, 226/226; differential vs 4b80565 on 212 inputs: only the new tests and post_apply1
  (registered before, now consumed). Found C90 (history file names) and C91 (force_edge_water
  -yes stored as a raw enum). Commit ea96bf6. Next: owner decisions for U2.
