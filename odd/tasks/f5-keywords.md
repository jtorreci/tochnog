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
- [x] T1 W1 switches and control (35): control_check_data, control_element_group,
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
  - [x] U2 aliases and run-level switches: control_reset_value_logarithmic_first (+ data.cc
    formula `log(b+c*x*d+e*x)` vs manual `ln((b+cx)/(d+ex))`: owner decision),
    global_element_dof_apply (`-only` new), message, print_data_name, control_check_data,
    control_zip, control_system_call (runs a user script: owner decision), processors_maximum
    (Professional default -yes vs GNU 1 thread: owner decision).
    control_system_call done `f6c4476` (owner: implement as the manual says).
    processors_maximum `1a4cb7d` (default -no until F6, owner); reset exponent/logarithmic_first,
    message, print_data_name, control_check_data, control_zip `f43142d`; global_element_dof_apply
    `03c8874`. Done.
  - [x] U3 post gates: control_post, control_post_apply (+ global post_apply, registered, unread),
    control_post_element_force. Done `ea96bf6`.
  - [x] U4 timestep loop: control_timestep_iterations_automatic_minimum_maximum_wished (GNU
    hardcodes 2/20/4, manual default 2 8 4), control_timestep_adjust_minimum_iterations,
    control_timestep_until_data / _maximum / _minimum.
    Done `9d0da65` (with U5).
  - [x] U5 element-group selection: control_element_group, element_group_apply,
    control_element_group_apply, mesh_element_group_apply, control_mesh_element_group_apply,
    processors_partition, print_precision.
    Done `9d0da65`.
  - [x] U6 input / database rewrite: mstring, mrange, control_input, groundflow_apply,
    linear_calculation_apply.
    Done `97c3e47` (delegated writer, context backstop; independent review PASS).
- [ ] T2 W2 output records (69): `element_*`, `post_*`, `print_*`, `control_print_*_node_type`,
  `materi_strain_*` / `materi_history_variable` initial dofs, `dof_*`, `global_*`,
  node_bounded_index (list in the 2026-11-01 file). Route: inline.
  Map (read-only explorer, 2026-11-02): 69 names; units:
  - [x] V1 pure aliases: control_print_{dof_line,gmsh,vtk}_node_type, element_contact_spring_direction,
    element_interface_intpnt_{strain,stress}_average, materi_history_variable, materi_strain_energy,
    post_node_rhside_ratio_dof_type. element_beam_force_moment moved to V4 (layout check).
  - [x] V2 global print defaults: print_{gmsh,vtk,tecplot}_calculation, print_gmsh_dummy,
    print_{gmsh,vtk}_node_type, print_vtk_coord, print_vtk_group, print_frd_freecad/prepomax
    (owner decision: manual default -yes for vtk/tecplot changes GNU output).
  - [x] V3 geometry/boundary records: element_geometry_present, print_element_geometry_present(_node_type),
    element_boundary (+ mesh_boundary). Inline. `aaf581c` (+ word -plus_displacement).
  - [ ] V4 element diagnostics: element_intpnt_{h,iso_coord,method,npoint}, element_truss_strain(_temperature),
    element_interface_intpnt_{direction,gap_status}, node_bounded_index, element_beam_direction_z,
    element_print_group_data_values.
  - [ ] V5 dof records/switches: dof_element_dof, dof_limit, global_node_dof_empty,
    global_element_dof_from_node_dof, element_intpnt_dof (decision: record vs GNU element storage).
  - [ ] V6 element_dof_initial_specific_number/value.
  - [ ] V7 materi initial dofs: materi_strain_plasti_{camclay,mohr_coul,tension,vonmises},
    materi_strain_total_{,compression_,shear_,tension_}kappa (decision: register only or fill).
  - [ ] V8 post_calcul modifiers: post_calcul_{absolute,label,limit,multiply,apparent_total,safety_default,static_pressure}.
  - [ ] V9 post data/point: post_count, post_data(_factor,_result), post_group_volume_summed, post_node_factor,
    post_node_rhside_ratio_method, post_point_{element_group,eps_iso,node_type},
    post_quadrilateral_element_group, post_strain_volume_initial.
  Professional examples only in validation_14_mesh.dat (output dump) and tochnog_abaqus.dat.
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
- 2026-11-02 owner on control_system_call: "dale, implementalo así" (as the manual: posix_spawn of
  ./tochnog_system_call.sh, no shell, warning when missing). Professional probes: one call per
  converged timestep (before its history print) + one at the end of each index with the record;
  stdout kept, stderr dropped, exit ignored, missing script silent, -3 rejected; the Professional
  only finds the script through the PATH (DIVERGENCES SYSCALL-PATH, class D). RED rc 1 on 9574b84;
  GNU call log identical. build_safe 376 OK / 0 FALLO, 226/226; audit_db_buffers rc 0. Also marked
  the U1/U3 keywords [x] in SEGUIMIENTO (missed before). Commit f6c4476. Next: rest of U2.
- 2026-11-02 owner on processors_maximum: "creo que mejor aquí converger con Pro. Para casos cortos
  casi no se va a notar en el flujo de trabajo general, y en casos largos es imperativo".
  Professional 02-08-2026 writes `processors_used 16` without a processors record (default -yes =
  every CPU). GNU threads were never exercised by the suite: build_safe with
  TOCHNOG_OPTIONS_PROCESSORS=4 gave 375 OK / 1 FALLO (force_point_sign1: disx 8 = 4 x 2,
  force_point_calculate ran inside the threaded parallel_new_dof_before, once per thread; moved to
  a serial call in top(); test force_point_threads1, RED rc 1 on the base binary). Known
  single-thread-only paths: group_user_umat (forced to 1), node_support_edge_normal_force
  accumulation (area.cc, warning). Corpus A/B 1 vs 4 threads on omen (mt_runs) before the default
  flips.
- 2026-11-02 owner on the threads measurement: "dale, esperá la suite y seguí" (keep the GNU default
  -no until F6). U2 closed: 1a4cb7d (processors_maximum, force_point serial; independent verifier
  PASS, 1-thread byte-identical), f43142d (reset exponent / logarithmic_first quotient = Professional,
  message, print_data_name, control_check_data, control_zip), 03c8874 (global_element_dof_apply,
  -only per manual; Pro aborts). build_safe 396 OK / 0 FALLO. Corpus A/B 1 vs 4 threads (omen
  mt_runs): pro_other 60/330 differ, 2 abort (interface_patch, patch1), ~13 differ > 1 % -> C92
  F6 prerequisite. Next: U4 (timestep loop).
- 2026-11-02 U4 + U5 (route inline, one commit 9d0da65: shared top.cc/database.cc). Professional
  probes: control_timestep_until_* stops the index after the step that leaves [min, max] and the
  next index runs (8 probes equal, missing maximum = error); automatic iterations default 2 8 4
  (Pro stops at 8; GNU had 20) but the step-size law differs (Pro x1.25 / x1.5625, /4; C94,
  class B); element_group_apply + control_element_group_apply = Professional; Professional reads
  mesh_element_group_apply without index and applies control_mesh_element_group_apply only at
  index 0 deleting the other elements (MESH-EG-APPLY, class D, GNU follows the manual), rejects
  control_element_group (CTRL-EG, class D); print_precision: manual example significant digits
  (PRINT-PREC, class C). RED: every new test rc 1 on the base binary (keywords unknown, 20
  iterations). build_safe 412 OK / 0 FALLO. U4 differential on omen (17 inputs with automatic
  stepping, mt vs u4): identical so far (10/17 done). Next: U6 (mrange, mstring, groundflow_apply,
  linear_calculation_apply, control_input; Pro examples tutorial_4, mesh_gid_* + 10.dat).
- 2026-11-02 U4 differential done (omen u4_runs, 17 automatic-stepping inputs, mt vs 9d0da65): 14
  identical; truss7 (index ends at t 0.069 after the step collapses, target still met), taylor3
  (differs, fails before and after), vajont_slope_fs (minimal-step stop before, > 900 s now) change
  with the 2 8 4 default; no Professional reference for them (parse / timeout). Recorded in C94.
  U6 delegated to one bounded writer (context backstop).
- 2026-11-02 U6 (route delegated: one bounded writer, context backstop; the NAS went down mid-unit,
  resumed after the owner restarted it: "he reiniciado el NAS. Podemos ocntinuar"). Professional
  probes in scratchpad/f5/u6/NOTES.txt: mrange counting = Professional (MRANGE-COMPILED class C);
  Professional rejects every mstring (MSTRING-PRO class D); groundflow_apply -no = Professional;
  linear_calculation_apply = manual list (Professional extras followed where the manual is silent;
  LINCALC-WOLFF, LINCALC-CONTACT class D); control_input = Professional (start of index, replace /
  add, errors; CTRL-INPUT-NEWELEM class C, C96). 22 new tests RED on f83f4ea; build_safe 431 OK /
  0 FALLO; differential 339 test-2026 inputs: only the new tests, timing and the C91 enum shift.
  Independent verifier (item 3, input_runtime refactor): PASS, 343-input differential identical
  (C91 shift and one concurrency flake), 13 runtime-file probes identical, ASan clean. Commit
  97c3e47. T1 closed. Next: T2 W2 output records (69).
- 2026-11-02 S5 refresh after W1: probe on the c135eeb objects: 34/320 now accepted, missing 286 =
  63 GiD + 223 (KEYWORD-GAP-REPORT update).
- 2026-11-02 T2 V1 (route inline): 9 aliases (db_number data part; initialization names
  normalised before they are stored). RED: alias_names_2026_1 rc 1 on c135eeb objects (unknown
  initialization word). build_safe 433 OK / 0 FALLO.
- 2026-11-02 owner on the V2 end-of-run prints (manual 6.998/6.999 default -yes for tecplot/vtk; the
  Professional 02-08-2026 writes no .vtk/.plt unless -yes is given): "vamos con pro, y lo dejamos
  documentado" -> GNU default -no (= Professional), class A divergence from the manual.
- 2026-11-02 T2 V2 (route delegated: one bounded writer, context backstop): 10 keywords + the
  control_print_frd_* records consumed; end-of-run prints default -no (owner); print_vtk_group adds
  the element_group cell field (manual = Professional; VTK-GROUP-NO class D for -no); 2D vtk/gmsh z
  garbage fixed. 22 tests; build_safe 455 OK / 0 FALLO. Independent verifier (item 3): PASS, 607
  inputs + 16 probes, only the element_group block and the z fix differ. Commit 05232b9. Open owner
  questions: C98 (print_vtk_coord meaning differs from the Professional), C99 (vtk/tecplot file
  names).

- 2026-11-02 owner on C98: "dale, seguí a Pro" -> print_vtk_coord / control_print_vtk_coord follow the
  Professional: POINTS always, VECTORS coord point field (stored coordinates) dropped by -no. Route
  inline. Commit 37d81d1; build_safe 455 OK / 0 FALLO. Next owner question: C99 (vtk/tecplot file names).
- 2026-11-02 owner on C99: "sigamos a Pro"; on the -separate_* formats (manual 6.339 <base><index>.vtk
  vs Professional <base>_000.vtk): "el de Pro". Route inline. Commit ecabcae: print_vt.cc / print_te.cc
  names = Professional (incl. -separate_index / -separate_sequential, .plt emptied per run), class A
  VTK-TECPLOT-NAMES; 6 new test-2026 tests (RED on e2e7c0c), test-2014 vtk_*1 on -separate_index;
  build_safe 461 OK / 0 FALLO. Next: T2 V3 (geometry/boundary records).
- 2026-11-02 owner on V3: "como el manual. Documentamos la divergencia" (mesh_boundary default -no;
  element_boundary only on elements with a free side). V3 committed `aaf581c`: area.cc mesh_boundary_set
  (renamed from contact_node_boundary_set), geometry.cc element_geometry_present_calculate, mesh.cc deletes
  element_boundary, database.cc registers -plus_displacement (rejected before in every record). Divergences
  MESH-BOUNDARY-DEFAULT (A), ELEMENT-BOUNDARY-ALL (D), EGP-PLUS-DISP (D); BACKLOG C40 (3) resolved.
  build_safe 475 OK / 0 FALLO, audit 0. Next: T2 V4 (element diagnostics).
- 2026-11-02 V4 Pro probes (scratchpad/f5/v4/in, omen ~/tochnog-runs/v4), Professional 02-08-2026:
  - element_intpnt_{h,iso_coord,method,npoint} and node_bounded_index are written in EVERY run (no switch;
    manual silent on the default -> follow Pro). quad4 2x2: method -gauss, npoint 4, h = per IP the nnol values
    (IP-major), iso_coord = per IP ndim values (+-0.57735); tria3: method 0, npoint 1, h 1/3 x3, iso 1/3 1/3
    (area coordinates); quad4 interface (group_interface -yes): -lobatto, 2 points, iso (-1,0),(1,0).
  - node_bounded_index: one slot per node_dof entry (npuknwn), the bounda_dof index of the record that
    bounds the dof, -nothing elsewhere (unbounded nodes: all -nothing).
  - element_truss_strain = total (length increase / length) 0.01; element_truss_strain_temperature =
    alpha dT 0.01 (alpha 1e-3, dT 10; the force is 0).
  - print_group_data -group_materi_elasti_young -group_materi_elasti_poisson -> Pro writes
    `element_print_group_data_value 1 1000 0.3` (SINGULAR; manual 6.446 `_values`). GNU print_group_data
    has length 1 (only one name parsed): must be widened.
  - element_interface_intpnt_direction / _gap_status: never written by Pro (also with an opened gap);
    Pro writes `element_interface_direction` (9 values per IP: normal, first and second tangential).
  - Beams: element_beam_direction_z overrides group_beam_direction_z per element (3D cantilever, Iyy 1 /
    Izz 4: tip vely 3.33e-3 with local z = global y, 8.33e-4 with z = global z); element_beam_force_moment 12
    values (first node minus). GNU lacks group_beam_direction_z(_reference_point) and
    element_beam_force_moment: beam family split to V4b. Pro 2D beam bounda_dof rot inputs rejected.
  - Not V4: GNU element_interface_intpnt_strain is the step increment (Pro total), BACKLOG line ~3438.
- 2026-11-02 V4a committed `ab1a014` (inline): element_truss_strain(_temperature), node_bounded_index,
  print_group_data list + element_print_group_data_values (= Pro; PGD-VALUE-NAME class D). build_safe 480 OK /
  1 FALLO (interface_tension_status1 time_calculation wall-clock flake, rerun identical), audit 0.
  Remaining V4: V4b element_intpnt_{h,iso_coord,method,npoint} (always written, Pro) +
  element_interface_intpnt_{direction,gap_status} (manual; Pro writes element_interface_direction only);
  V4c beam family (group_beam_direction_z(_reference_point), element_beam_direction_z,
  element_beam_force_moment).
- 2026-11-02 V4b committed `71adae7` (delegated writer, context backstop): element_intpnt_{h,iso_coord,method,npoint}
  via element_intpnt_set (polynom.cc, step_start) + pol() thread-local method/iso outputs; Pro-equal for bar2/3,
  quad4/9, tria3/6, hex8, tet4/10, 5 interface shapes; new class B rows INTPNT-H-LENGTH, TET4-MAX-ORDER (tet4
  -maximal points 3/4 swapped vs Pro: reorder needs an owner decision). build_safe 515 OK / 0 FALLO, audit 0;
  parent spot check element_intpnt_quad4_1 = Pro. Independent verifier (Risk item 3: pol(), bounda(),
  mesh_boundary_set shared code) deferred to one run at the end of T2 covering V3, V4a, V4b.
- 2026-11-02 V4c committed `e17a202` (delegated writer): beam_local (group/element_beam_direction_z, reference
  point), element_beam_force_moment (12 local values; planar beam too). build_safe 580 OK / 0 FALLO, audit 0;
  15 existing beam inputs unchanged apart from the new record; parent spot check element_beam_direction_z1.
  Open owner decisions (class B): C37 (3D beam without record stays planar vs Pro full 3D), BEAM-MEMORY-DEFAULT
  (GNU -updated vs Pro -total_linear: 0.2399880 vs 0.24), BEAM-DIRECTION-LENGTH (2 vs 9 values), and
  TET4-MAX-ORDER (owner asked about side effects; recommendation: reorder). Pending: V4 interface
  element_interface_intpnt_{direction,gap_status}.
