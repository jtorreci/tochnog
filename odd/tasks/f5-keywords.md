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
- S6 New elements (owner, 2026-11-09, U9): add prism12, prism18 and a continuum quad6 so that
  control_mesh_extrude_element and control_mesh_convert_quad9_quad6 are implemented, verbatim: "A mi me parece
  que es interesante añadir los elementos". VTK output is included for the new elements and prism15, verbatim:
  "si, la incluimos. Si vtk no lo soporta directamente convertimos en la escritura del vtk a un tipo soportado".

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
- [x] T2 W2 output records (69): `element_*`, `post_*`, `print_*`, `control_print_*_node_type`,
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
  - [x] V4 element diagnostics: element_intpnt_{h,iso_coord,method,npoint}, element_truss_strain(_temperature),
    element_interface_intpnt_{direction,gap_status}, node_bounded_index, element_beam_direction_z,
    element_print_group_data_values.
  - [x] V4d beam owner decisions (2026-11-02): native 3D beam for every 3D beam (C37), default
    group_beam_memory -total_linear (BEAM-MEMORY-DEFAULT), element_beam_direction 9 values
    (BEAM-DIRECTION-LENGTH). Route: inline (no writer reason). Commit `f99793a`.
  - [x] V5 dof records/switches: dof_element_dof, dof_limit, global_node_dof_empty,
    global_element_dof_from_node_dof, element_intpnt_dof. Owner decision 2026-11-02: option 2, rename
    the GNU per-integration-point storage to element_intpnt_dof and make element_dof the per-element
    nodal values (manual 6.421 / 6.437 = Professional). Persisted format change (Risk item 1):
    own commit + independent verifier. Route: mapping explorer first, then decide writer/inline.
  - [x] V6 element_dof_initial_specific_number/value. Route: inline (no writer reason). Commit `df07925`.
  - [x] V7 materi initial dofs: materi_strain_plasti_{camclay,mohr_coul,tension,vonmises},
    materi_strain_total_{,compression_,shear_,tension_}kappa. Owner 2026-11-03: follow the manual (fill),
    document the Professional divergence. Route: inline (no writer reason). Commit `0fa067b`.
  - [x] V7b C107: older per-model plastic strain dofs (hardsoil/cap/compression/diprisco/druckprag) gated
    by model with the epp recovery (0 under another model). Owner 2026-11-04. Route: inline. Commit `1a44872`.
  - [x] V7c C108 labels: Professional dof labels for the per-model plastic strains and the kappas. Owner 2026-11-04.
    Route: inline. Commit `c854d82`.
  - [x] V7d C108 kappa measures: compression sign (Professional negative), history per manual; doc corrections.
    Route: inline. Commit `94c6161`.
  - [x] V8 post_calcul modifiers: post_calcul_{absolute,label,limit,multiply,apparent_total,safety_default,static_pressure}.
    Route: inline. Commit 42dac7d. C109 decided: keep the GNU manual names (owner 2026-11-04, "mantenemos nombres GNU-manual").
  - [x] V9 post data/point: post_count, post_data(_factor,_result), post_group_volume_summed, post_node_factor,
    post_node_rhside_ratio_method, post_point_{element_group,eps_iso,node_type},
    post_quadrilateral_element_group, post_strain_volume_initial.
  Professional examples only in validation_14_mesh.dat (output dump) and tochnog_abaqus.dat.
    Route: inline (explorer for the post.cc map). Commit a1c7d7d. C110 open for the owner; C111 decided (manual, class A).
- [x] T3 W3 materials (27): `group_materi_*` and `materi_plasti_*` (umat, stress_null, hypo
  niemunis/wolffersdorff variants, borja_tamagnini, elasti_c, hyper reduced_polynomial /
  simo_taylor, visco_power, residual_stiffness, mpc, diprisco_density, history_variable_user,
  phimob, maximum_iterations, hypo_substepping). Route: per family.
  Probe 2026-11-06: 29 rejected (24 group_materi_* + materi_history_variable,
  _strain_energy, _plasti_hypo_substepping, _plasti_maximum_iterations, _plasti_phimob);
  explorer map in Log L-T3-map. Family order A, B, C1-C3, D1-D5, E.
  - [x] A materi_history_variable / materi_strain_energy: already normalised in input()
    (initialization part, 60b30d0); the db_number probe cannot see it. No change.
  - [x] B aliases group_materi_hyper_reduced_polynomial, _volumetric_simo_taylor,
    group_materi_umat(_parameters) + global materi_plasti_maximum_iterations. Route: inline.
    Commit c095d18. C112 resolved (default 20000, owner), commit 6a78839.
  - [x] C1 hypo_substepping (+ unused control record), hypo_minimum_void_ratio. Route: inline.
    Commit 6745219.
  - [x] C2 plasti_mpc + mpc_factor. Route: inline. Commit a139f16. C114 (wall cohesion) found.
  - [x] C114 wall reduction of c and the direct-path factor (manual + Pro agree). Route: inline.
    Commit 397a3cb.
  - [x] C3 residual_stiffness, phimob. phimob: route inline, commit 5e42bae. residual_stiffness: route
    inline, commit a806c1c (owner 2026-11-06: manual followed, DIVERGENCES PLASTI-RESIDUAL-STIFFNESS).
  - [x] D1 elasti_c + c_direction (route inline, commit 6b8327d).
  - [x] D2 visco_power_name / _value (route inline, commit 2ef6bdc).
  - [x] D2x visco_exponential_name / _values / _limit (route inline, commit 8f8df2a; owner 2026-11-07: manual
    followed, DIVERGENCES VISCO-EXPONENTIAL-NAME; C117 found).
  - [x] D3 wolffersdorff_niemunis (route inline, commit 46a0377). Prerequisite C118 (void_ratio_linear
    alias) resolved, route inline, commit 932c5f9.
  - [x] D4 borja_tamagnini (route inline, commit df32477; owner 2026-11-07: paper over Pro, DIVERGENCES
    BORJA-TAMAGNINI).
  - [x] D5 stress_null family (route inline, commit 1fc60e6).
  - [x] E1 user hooks young_user, history_variable_user (+ _parameters) as stubs with the Professional errors,
    umat_pardiso_decompose without effect (class A). Route inline, commit 536998f.
  - [x] E2 diprisco_density (route inline, commit 6059c82; owner 2026-11-07: manual over Pro, linear
    interpolation with the current void ratio, DIVERGENCES DIPRISCO-DENSITY; C119 found).
  - [x] E3 niemunis_visco (+ _ocr, control ..._ocr_apply): route inline (explorer mapped the thesis), commit
    8f45881; owner 2026-11-07: b2 < 0 rejected per manual (DIVERGENCES NIEMUNIS-VISCO-B2); C120 found.
- [x] T4 W5 geometry, mesh tools, mpc, dependency (69). Route: per family. Closed 2026-11-13 (verifier, 16f3721).
  Probe 2026-11-08 (clean build of cdcf986): exactly 69 W5 names rejected. Map by a read-only explorer
  (mapping trigger: 69 names across ~12 files), summary in Log L-T4-map. Units:
  - [x] U1 aliases control_mesh_delete_geometry_element_group, geometry_triangle_eps_iso,
    mesh_interface_triangle_coordinates. Route: inline. Commit 87436f7.
  - [x] U2 geometry entity options (12): geometry_{line,quadrilateral}_eps_iso, geometry_boundary, geometry_mpc,
    geometry_exclude, geometry_element_geometry(_method), geometry_circle_part, geometry_cylinder_part(_start_vector),
    geometry_tetrahedral, geometry_hexahedral. Route: inline (one file family, no writer reason). Commit 073d4f3.
    geometry_method moved to U4 (element selections).
  - [x] U3 dependency modifiers (4): dependency_{method,type,number,geometry} (+ element group -all of the manual example;
    print_group_data through get_group_data). Route: inline (no writer reason). Commit 3b66fc9. C122 found; decided 2026-11-08 class D, GNU kept.
  - [x] U4 delete/refine flags (7): geometry_method, control_mesh_delete_geometry_{direct,method,projection_type},
    control_mesh_refine_locally_{minimal_size,not_method,only_method}. Route: delegated (context backstop).
    Commit bbbc741. C123 (projection default, owner) and C124 found.
  - [x] U5 mesh operations (15 names; the "16" counted control_mesh_extrude_element, which is U9). Map: explorer
    L-T4-U5-map. Sub-units, all inline (shared database.cc / headers / top.cc / build_safe, so no disjoint
    parallel surfaces):
    - [x] U5a extrude (d040c57): control_mesh_extrude layer convention fixed to manual 6.179 (z0..zn layers, each split by
      _n; GNU used the untouched z=0 nodes as an extra bottom layer and ignored _n for linear elements), plus
      control_mesh_extrude_direction, _element_group_new, _contact_spring_element_group(_new).
    - [x] U5b split / convert (140e7a0): control_mesh_split_element_from, _split_element_to (-tet4 / -prism6),
      control_mesh_convert_tria6_tria3. control_mesh_convert_quad9_quad6 moved to U9 (no GNU quad6 continuum
      element).
    - [x] U5c (394ac2a) control_mesh_renumber_element_geometry_offset / _group_offset,
      control_mesh_duplicate_element_group, control_mesh_remove_geometry. control_mesh_remove_keep_geometry moved to U9
      (only meaningful with control_mesh_remove_really -no: the GNU deletes removed elements).
    - [x] U5d (ea220df) control_mesh_generate_truss_beam_separate, mesh_correct.
  - [x] U6 (5db4dd0) mpc_element_group options (6): control_mpc_element_group, mpc_element_group_{always,closest,
    coord_geometry,keep,eps_iso}. Route: delegated (bounded writer).
  - [x] U7 (663cbcb) control_mesh_truss_distribute_mpc family (10 records, 6.245-6.254) and mpc_geometry
    _method/_tolerance with the parents mpc_geometry/_switch/_dof (6.866-6.870). Owner: implement the full feature.
    Route: delegated (bounded writer).
  - [x] U8 (dd2ee02) geometry_moving family (6) + control_geometry_moving: new module geometry_moving.cc.
    Route: delegated (bounded writer).
  - [x] U9 owner decisions: all taken (2026-11-11).
  - [x] U9f (d13d4e2, 728e997) (S2-S5) (C126) control_mesh_remove -method1 geometric containment (manual 6.229 "completely located inside")
    and -method3 (all nodes have node_mpc). Owner: implement (2026-11-13). Route: delegated (bounded writer, same as
    U9e; context backstop). RED: 22 of 40 new checks fail on 31c00c6; GREEN suite 1231 OK, 0 FALLO. Integration points
    (owner 2026-11-13, 728e997): RED gap2 / lcorner1 fail on d13d4e2, GREEN suite 1235 OK, 0 FALLO; class D
    REMOVE-METHOD1-INTPNT. REMOVE-METHOD1-CURVED class D confirmed 2026-11-13; tolerance values kept (tet / prism not measured).
  - [x] U9e (5f08f81, 9c08983) (S2-S5) control_mesh_remove_really (6.232), _really_activate_all (6.233), _really_activate_factor
    (6.234), control_mesh_remove_keep_geometry (6.231). Owner: implement (2026-11-11). Route: delegated (bounded writer;
    context backstop of the parent session). RED: base 896978f "I don't know what to do with : control_mesh_remove_really",
    32 of 34 new checks fail; GREEN suite 1182 OK, 0 FALLO. Nodal start (owner 2026-11-11, 9c08983): RED 3 of 43 U9e
    checks fail on ac4bc9b, GREEN suite 1191 OK, 0 FALLO; class D REMOVE-REACTIVATION-START with the recipe (factor 0 +
    _activate_all -yes = Professional, 5e-11). Owner questions open: REMOVE-INACTIVE-RECORDS (B),
    REMOVE-NO-TIMESTEP / -ACTIVATE-ALL / -PRESSURE-GRAVITY (D), -method1 geometric containment (C126).
  - [x] U9d (b3f4746) (S2-S5) control_mesh_map (manual 6.210, Pro test/other/map1): quadratic <-> linear mesh map as a
    solver first guess. Owner: implement (2026-11-11). Route: delegated (bounded writer; context backstop of the parent
    session). RED: base stops "I do not know : control_mesh_map"; GREEN suite 1148 OK, 0 FALLO. Owner decisions 2026-11-11: C134 = Professional,
    MESH-MAP-TET10 / -PRISM15 class D, MESH-MAP-MID-NODES class B kept. C134 decided 2026-11-11 (= Professional).
  - [x] U9a (3c93898, a1185ed, 505564d, b006671) (S3, S4, S6) Professional probes for quad6, continuum quad6,
    control_mesh_convert_quad9_quad6, VTK output of quad6 and prism15. Route: delegated (bounded writer; context
    backstop). Suite 1086 OK, 0 FALLO. QUAD6-3D class D confirmed by the owner (2026-11-10).
  - [x] U9b (0b35e26, c0f8038, b9599d0, af188d0) (S2-S6) prism18, control_mesh_extrude_element with -tria6 (manual
    6.181 default -prism18), VTK cell 32. Route: delegated (bounded writer; context backstop, the parent session was
    already compacted once). Suite 1102 OK, 0 FALLO. PRISM18-MINIMAL class A confirmed by the owner (C132, 21 points).
  - [x] U9c (fd32b82, e6fa737, dd87da9, 1c258be) (S2-S6) prism12 continuum, -prism12 option of the extrude, VTK cell
    31. Interface prism12 is out of scope (refused, class B). Route: delegated (bounded writer; context backstop, as
    U9a and U9b; 6104401 quad9 fix). Suite 1119 OK, 0 FALLO. Owner confirmed PRISM12-DEFAULT / PRISM12-POST-POINT class D (2026-11-11).
- [x] T5 W6 other physics (25): baseline correction, force_edge_diagram, force_gravity_geometry,
  force_volume_element_group_0, groundflow_phreatic_*, group_beam_*, group_contact_spring_*,
  group_dof_initial*, group_time_fill, group_truss_elasti_elongation_force_diagram,
  slide_damping, slide_user, strain_volume_*. Route: per family.
  Map 2026-11-13 by a read-only explorer (mapping trigger: 24 names across ~12 files), summary in Log L-T5-map.
  Units:
  - [x] U1 aliases group_contact_spring_plasti_cohesion, group_contact_spring_plasti_friction_automatic,
    groundflow_phreatic_level_n. Route: inline. RED: base "I do not know : group_contact_spring_plasti_cohesion";
    GREEN alias_names_2026_4 (Professional stores the same values).
  - [x] U1b groundflow_phreatic_bounda -yes/-no (6.570). Route: inline. Owner decision L-T5-bounda.
    RED: base "Error in data part" on the name. GREEN groundflow_phreatic_bounda1/2; default -yes kept
    (absent = -yes measured on the Professional, and the GNU free-surface block already matched it).
    Divergence PHREATIC-BOUNDA-NO-POST.
  - [x] force_volume_element_group_0 discarded (PDF header artifact of 6.504; owner L-T5-fveg0). Gap report: 57 pending.
  - [x] U3a groundflow_phreatic_only (6.581). Route: inline. Owner L-T5-only. RED: base "I do not know :
    groundflow_phreatic_only". GREEN groundflow_phreatic_only1/2/_error1 (Professional .dbs: group_type 2 -materi,
    group_groundflow_permeability 2 deleted; total pressures equal).
  - [x] U3b inclined single level + groundflow_phreatic_project (6.582). Route: inline. Owner L-T5-project.
    RED: base inclined1 misses the targets (no seepage), "I don't know ... groundflow_phreatic_project".
    GREEN groundflow_phreatic_inclined1, _project1/2, _project_error1 (Professional totals to 1e-6).
    Divergence PHREATIC-PROJECT-3D.
  - [x] U3c inclined groundflow_phreatic_level_multiple (same fix per record; _static -yes untouched). Route: inline.
    Owner L-T5-multiple. RED: groundflow_phreatic_multiple_inclined1 misses the targets (dof 0, no seepage); GREEN
    with the Professional totals to 1e-6.
  - [x] Owner decisions: force_edge_diagram (resolved by measurement, L-T5-diagram),
    group_contact_spring_direction_automatic_planes meaning (resolved by measurement, L-T5-planes),
    group_beam_force_moment_plasti (resolved, L-T5-beamplasti-owner), strain_volume kernel (resolved, L-T5-strainvol-owner), baseline correction scope (resolved, L-T5-baseline-owner).
  - [x] U2 slide_damping (viscous, divergence SLIDE-DAMPING-VISCOUS), control_slide_damping_apply, slide_user
    (user.cc stub, Professional message). Route: inline. Owner L-T5-slide. RED: base "I do not know : slide_damping"
    / "slide_user". GREEN slide_damping1 (c*v exact), slide_damping2 (Professional values), slide_user_error1.
  - [x] U4 group_dof_initial, _specific_number, _specific_value. Route: inline. Owner L-T5-dofinit (option B,
    element records first; divergence DOF-INITIAL-GROUP-PRECEDENCE). RED: base "I don't know what to do with :
    group_dof_initial". GREEN group_dof_initial1..7 + _error1 (Professional values except _4).
  - [x] U5 force_edge_diagram (per edge node factor from the vertical displacement of the iterate, 0 outside),
    force_gravity_geometry (elements with all nodes inside; materi and truss). Route: inline. Owner "si" (measure
    first), L-T5-diagram, L-T5-gravgeo; no divergence. GREEN force_edge_diagram1/2/_error1,
    force_gravity_geometry1/2/3 (Professional values).
  - [x] U6 group_contact_spring_direction_automatic_planes (a -no switch removes that component of the automatic
    direction, nothing left = error) plus three measured contact spring fixes (no friction record = unlimited
    friction, 2D in-plane tangent, signed matrix); group_truss_elasti_elongation_force_diagram (force from the
    diagram of the total strain, 0 outside, young and initial force ignored). Route: inline. Owner "si" (measure
    first), L-T5-planes, L-T5-trussdiag. RED: base "Error in data part" (both keywords), sideways conspr disx 1.5e12.
    GREEN contact_spring_tangential1, _inclined1, _planes1/2, _planes_error1, truss_force_diagram1..5, _error1.
  - [x] U7 group_time_fill (alive birth_empty..death; density and total pressure times f = (t - be)/(bf - be) at
    the step end; no plasticity while f < 1). Route: inline. Owner "si" (measure first), L-T5-timefill. RED: base
    "Error in data part"; plasticity skip off: group_time_fill2 disy -0.01153; pressure unscaled: -0.0025.
    GREEN group_time_fill1/2/3 (Professional values).
  - [x] U8 strain_volume_* (6.1062-6.1066) + post_strain_volume_absolute / _relative. Route: inline (one behavior,
    no writer reason). Owner L-T5-strainvol-owner. RED: base "Error in data part" (relative_time, element_group,
    geometry, list) and no post results. GREEN: 10 tests, Professional values to 1e-8; divergences
    STRAIN-VOLUME-ELEMENT-LIST, STRAIN-VOLUME-BOTH-TABLES. Commit a8ed5db.
  - [x] U9 bounda_baseline_correction (6.19) + _parameters (6.20). Route: inline (one behavior, no writer reason).
    Owner L-T5-baseline-owner (manual over the Professional). RED: base "I don't know what to do with :
    bounda_baseline_correction". GREEN bounda_baseline_correction1/2/3 + _error1; divergence
    BASELINE-CORRECTION-SIGNAL. Commit ea50f33.
  - [x] U10 group_beam_force_moment_plasti (6.589): plastic hinges at the beam element ends (return mapping of the
    end moments per bending plane, shear and torsion limits, elastic unloading). Route: inline (one behavior, no
    writer reason). Owner L-T5-beamplasti-owner (real model, not the Professional clipping). RED: base "I do not
    know : group_beam_force_moment_plasti". GREEN group_beam_force_moment_plasti1..5, _2d1, _error1 (analytic);
    divergence BEAM-PLASTIC-HINGE. Commit 31affab; suite 1301 OK / 0 FALLO, 226/226. RDD assess high (build_safe.sh
    shell signal), STATUS stop rdd_disabled: no review.
  - [x] End of T5: independent verifier (risk item 3: shared materi / bounda / conspr / beam code), baseline 9777b2d.
    1 blocker (signed zero in groundfl.cc phreatic bounds, U3b/U3c) fixed; all units PASS (L-T5-verify).
- [x] T6 W7 safety_slip (32): self-contained slip-surface safety-factor feature, planned
  separately. Route: own feature document odd/tasks/f5-safety-slip.md (2D first, 3D second, owner 2026-11-22).
  Commits 0964a37 (U1), 900f94a (U2+U3), 38e181e (U4+U5), 5aa6194 (3D measurements), 7aec39b (U6); suite
  1335 OK / 0 FALLO; safety 31/31, control_safety_slip and control_print_gid_safety_slip_critical done.

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
- 2026-11-02 Owner on the V4c beam decisions (verbatim): "1. Me gustaría tener vigas 3D nativas, 2. -updated
  se vendía en su momento como una feature de Tochnog; de hecho recuerdo que me llamó mucho la atención y tuve
  que corregir mis inputs y usar total-linear para dar los resultados clásicos de resistencia de materiales;
  quizá mejor -total_linear, pero documentando bien la opción. 3. Harán falta los 9 valores si queremos
  realmente una viga 3D nativa, no?" -> V4d: every 3D beam through beam_local with the Professional default
  local z (group_beam_plane, GNU-only, then gives the local z as the plane normal); default memory
  -total_linear with -updated documented; element_beam_direction = 9 values (local x, y, z rows, Pro layout
  probe dfltx/obl). 2D: Pro 2D beam singular on the probes, 9 values written by analogy (z = global z).
  TET4-MAX-ORDER answer still pending.
- 2026-11-02 V4d committed `f99793a` (inline): beam_3d() sends every 3D beam to beam_local() (record 0 = default
  local z; group_beam_plane normal if given, error if the beam is along it), memory default -total_linear,
  ELEMENT_BEAM_DIRECTION 9 values (beam_local R rows, beam_2d c s 0 -s c 0 0 0 1). RED on 3ba1901, GREEN;
  build_safe 634 OK / 0 FALLO, audit 0. trubea2/trubea3 = Professional (C37 resolved in BACKLOG). A/B: inputs
  without group_beam_memory move to the exact classical answer; explicit memory unchanged. DIVERGENCES: C37,
  BEAM-DIRECTION-LENGTH removed (= Pro), BEAM-MEMORY-DEFAULT -> class C BEAM-MEMORY-UPDATED, new class C
  GROUP-BEAM-PLANE. Risk item 3 (shared beam path, existing results change by owner decision): covered by
  the end-of-T2 independent verifier. TET4-MAX-ORDER still pending.
- 2026-11-02 Owner: "reordenemos, pues." -> TET4-MAX-ORDER committed `3dc43db` (inline): pol() tet4 -maximal points
  3 and 4 in the Professional order; element_intpnt_3d1 expectation updated; 634 OK / 0 FALLO; A/B on the two
  tet4 -maximal inputs: only element_intpnt_* and element_dof (points 3/4 exchanged) differ. Divergence row removed.
- 2026-11-02 Owner: "sí, como el manual, y documentamos la divergencia" -> V4 interface records committed `4aeb107`
  (inline): element_interface_intpnt_direction (9 per interface intpnt, normal / first / second tangential)
  and element_interface_intpnt_gap_status (strain_eff <= gap closed). Pro probes v4i d2/d3/g2 (Pro writes
  neither; element_interface_direction same vectors in 2D, hex8 normal reversed on half the points).
  RED on e94f261, GREEN; 643 OK / 0 FALLO, audit 0. Divergence INTERFACE-INTPNT-DIRECTION-GAP (class A).
  V4 complete. Next: V5 (decision element_intpnt_dof record vs GNU element storage).


- 2026-11-02 V5 element_intpnt_dof probe (Pro v4/r_intq quad4 2x2): Pro writes element_intpnt_dof (4 IP x 10)
  and element_dof (4 nodes x 10, the per-element nodal fit of miscel.cc); GNU element_dof = the IP values
  (270 = 27 x 10 padded), equal to the Pro element_intpnt_dof except 2 displacement slots (max 0.006).
  Owner (verbatim): "prefiero la 2, es como convergemos con manual y pro. Hacer el uno por retrocompatibilidad
  solo sería compatibilidad con el dbs, que es un archivo de salida, no de entrada". Note given: input.cc:1965
  also reads element_dof from a data part (a .dbs fed back as input keeps the old layout).
- 2026-11-03 V5 element_intpnt_dof committed `6467dea` (delegated writer, context backstop; mapping explorer before):
  ELEMENT_INTPNT_DOF = former IP storage (npoint*nuknwn, not padded; zero records of beam/spring/truss elements
  deleted as the Professional), ELEMENT_DOF = per-element nodal (miscel.cc element_dof_nodal_set, serial after
  the step-end promotion: principal dofs = NODE_DOF, IP states = recovery fit kind / fallback h*V mean).
  = Professional within 4e-11 on intq, two, p1, q9 except disx/disy (ELEMENT-DOF-DISPLACEMENT, class B).
  RED on bb69dd0, GREEN; 654 OK / 0 FALLO, audit 0; A/B on 739 suite inputs: rc identical, 543/545 .dbs
  identical outside the two records (2 gain a memory-report line), 12679 old element_dof = new element_intpnt_dof.
  Docs: 26 manual pages renamed ELEMENT_DOF -> ELEMENT_INTPNT_DOF (code token) + 4 IP-storage sentences.
  Residuals: interface/empty elements keep the 27-point zero element_intpnt_dof; an old .dbs fed back as input
  loses the IP state (element_dof_initialised 1, documented). Pre-existing nonloc.cc PUT outside
  if(ielem_nonloc) noted. Risk items 1 and 3 -> independent verifier next.
- 2026-11-03 V5 independent verifier (fresh mktemp copies of bb69dd0 / 6467dea, ASan on 16 inputs): PASS, no
  blockers. Probes = Professional within 3e-12 (dis slots excluded); 651 corpus inputs: rc identical, all other
  .dbs records identical (bicg_csr_threads1 1e-10 thread nondeterminism, reproduced on the baseline); 13180 old
  element_dof = new element_intpnt_dof; continuation from a new .dbs = continuous run. Advisories: old-format
  .dbs input loses the IP state (documented residual); interface zero record length = continuum rule (doc
  corrected); structural zero records removed only at a step end; -no keeps the zero record; pre-existing
  nonloc.cc PUT outside if(ielem_nonloc).
- 2026-11-03 V5b committed `cde109f` (delegated writer, context backstop): dof_element_dof (= Pro except work/f/void,
  C103), dof_limit (pairs per dof_label entry, node_dof clamped, = Pro; C104 prescribed-pass residual),
  global_node_dof_empty (default -yes, displacement kept as the Professional; mesh_element_group_apply groups
  not counted as empty), global_element_dof_from_node_dof (switch; default kept -yes). 688 OK / 0 FALLO,
  audit 0; A/B 748 inputs: only the new tests and global_node_dof_empty1/3 change. Owner decisions open:
  C100 default of global_element_dof_from_node_dof (manual + Pro -no, GNU always -yes, 214 suite inputs with
  node_dof), C101 global_node_dof_empty displacement (manual zero, Pro keeps). Pre-existing C102, C105 noted.
  Risk item 3 -> covered by the end-of-T2 independent verifier.

- 2026-11-03 Owner on C100: "como manual y Pro". Default of global_element_dof_from_node_dof flipped to -no.
  First flip broke hypo_masin_history_slots1 (initial state by control_reset_dof only; the Professional runs it
  — minipc Pro run): data.cc control_reset_dof now writes the integration points under -no. A node-zeroing
  variant was rejected (it wiped the input displacement, global_node_dof_empty1). Owner: "si" to writing -yes
  explicitly in the 20 node_dof-initialised suite inputs (+ convected_stress_recovery1); the Professional parser
  rejects them for other keywords. New test global_element_dof_from_node_dof3 (default). Suite 692 OK / 0 FALLO,
  226/226 runs, audit 0. Pending: A/B on the 748 inputs on omen (offline). Risk item 3 (shared reset code in
  data.cc; default change) -> end-of-T2 independent verifier.
- 2026-11-03 A/B omen (pre-C100 binary vs 3482d4b, pre-C100 inputs): 21 rc changes (the -yes inputs) and 4
  more silent node_dof-initialised inputs (mcap1, mcap1_comb, mstrain_cap, groundflow_total_pressure_tension)
  -> -yes. Reset at the integration points was an element average; the Professional evaluates the reset at the
  point coordinates (minipc Pro runs) -> data.cc reset_value_at() lambda + pol(); = Pro to 1e-12. New checks
  (exponent1, logarithmic_first1 element_intpnt_dof; global_node_dof_empty3 node 2 = Pro). C106 found
  (pre-existing shared diagram buffer). Owner on C101: "Me gusta más la idea de Pro y actual de GNU. Si los
  desplazamientos no están en cero y activamos esos elementos, corremos el riesgo, en una configuración
  -updated de que los nodos estén "dentro" de otros elementos" -> displacement kept, class A. Suite 694 OK /
  0 FALLO, 226/226, audit 0; A/B with current inputs: rc identical, 538 identical, 11 expected diffs.
- 2026-11-03 owner: "perfecto. Sigue". T2 V6 (inline): Professional probes (phase1 bar + quad4) show the
  field value_0 + grad . x evaluated at the nodes, birth step only, element_dof_initial taking precedence;
  constant value_0 alone NaN in the Professional (manual followed, DIVERGENCES DOF-INITIAL-SPECIFIC-CONSTANT).
  general.cc node field + top.cc marker; tests element_dof_initial_specific1..5 + _error1. RED on fe664a3.
- 2026-11-03 owner on V7 ("Pro accepts them as no-op and rejects camclay; manual says added to node_dof"):
  "yo creo que seguir el manual y documentar la divergencia con Pro". T2 V7 (inline):
  - Per-model plastic strain: filled only for groups with the model, with the epp recovery
    (materi.cc strain recovery block; general.cc weight list).
  - Total strain kappas: integration-point history maxima of new_ept (elem.cc masks), with
    materi_strain_total required.
  - Tests: materi_strain_plasti_model1..4 (equal to epp / 0); materi_strain_total_kappa1 (analytic);
    _error1.
  - Found BACKLOG C107: the older per-model dofs ignore the model and the epp recovery.
- 2026-11-04 owner on C107: "lo que hace Pro me parece más lógico. Si un elemento no tiene un modelo no
  tiene sentido rellenar sus dof, salvo que esos dofs también hagan falta para los modelos que sí que tiene
  ese elemento. [...]" then "sí, implementalo así. Los llevamos a cero si no tienen ese modelo".
  T2 V7b (inline): the 5 older dofs moved to the materi.cc per-model block (9 models) and the general.cc
  recovery weight list; old h_lump block removed. No model reads these dofs (output only).
  RED materi_strain_plasti_model5 (old binary: eppcap* = epphs* = 0.913406 under druckprag).
  Professional on model5: identical values (labels eppdr*/eppca*/eppha*). Suite 720 OK / 0 FALLO, 226/226.
  Found BACKLOG C108: the V7 "Professional adds no dof" reading was wrong (v7pro dbs hold eppvo*, kapto,
  strtocokap, ...); per-model values equal GNU, labels and kappa measures differ. V7 docs to correct.
- 2026-11-04 owner on C108 labels: "si, pásalas". Renamed 52 basenames in database.cc, mstrain_* targets, manuals.
  Pro analytic kappa probe (materi_strain_total_kappa1 input): kapto = history max sqrt(0.5 e:e) (= GNU);
  strtocokap -0.01 then -0.005 after unloading (current, signed); strtoshkap current (0.00721 -> 0.00361);
  strtotekap 0 alone, copies the shear value when both are initialised. Manual says "as occurred in history":
  S2 -> history per manual, compression sign from the Professional (manual silent on sign).
  V7d: compression kappa = signed history minimum min(old, min(eps_min,0)); RED = previous expectation 0.01
  (passed with c854d82), GREEN -0.01. DIVERGENCES STRAIN-KAPPA-HISTORY added, STRAIN-PLASTI-MODEL-DOFS
  narrowed to camclay; V7 docs corrected. Suite 720 OK / 0 FALLO, 226/226. C108 resolved.
- L-V8 (2026-11-04) owner: "pues venga, a por V8". Professional 02-08-2026 probes (scratch v8pro, v8pro2):
  multiply/limit per dof/operat item (3 or 5 factors for 3 items identical; 2 factors / 1 pair rejected
  "Error detected for data item : post_calcul_multiply/_limit"); order multiply -> limit -> absolute
  (q1 limit [0.1,0.2] + abs on -0.333 -> 0.1; q4 factor -2 + abs -> 0.667). apparent_total ground17:
  1.333333e7 / 1/3. safety_default ground14: 99 above the level. static_pressure -7 ground14: st -7,
  dy = total + 7, piping p_dyn too, lifting unchanged. post_calcul_label Pro `-a_ept -va0_ept`, GNU names
  `aept` (manual example) -> DIVERGENCES POST-CALCUL-LABEL-NAMES, BACKLOG C109 (owner decision).
  RED: old binary "Error in data part" on the 7 inputs. GREEN: node_dof_calcul = Pro (prival order K12).
- L-V9 (2026-11-04) V9 after C109 ("mantenemos nombres GNU-manual"). Explorer map: post() in post.cc, locator
  parallel_post_point (NODE_START_REFINED, point_el eps 1.e-3), post_node_rhside_fixed_free (max |.|, threshold
  1.e-6), strain_volume_* registered but not computed. Professional probes v9pro/v9t (3 quad4, 2 groups): count
  3 8 1, data 6, node factor x3, group volume 2 1, strain volume initial 1, point group, quadrilateral group
  0.0880956 equal. Pro ignores post_point_eps_iso and post_point_node_type -node -> DIVERGENCES POST-POINT-SEARCH,
  BACKLOG C111 (owner). rhside fixed/free maxima vs manual sums -> BACKLOG C110 (owner). RED: old binary rejects
  the records (r9 rc=1). GREEN: 5 new tests, targets met.
- L-C111 (2026-11-05, owner, verbatim): "manual+GNU, creo que es lo mejor" -> C111 decided: GNU keeps the manual
  behaviour for post_point_eps_iso / post_point_node_type; DIVERGENCES POST-POINT-SEARCH class A. Next: end-of-T2
  independent verifier.
- 2026-11-05 end-of-T2 independent verifier (fresh mktemp copies of 3227efe / 35daf7e, 612-input corpus
  differential, suite 744 OK, audit 0): PASS on V3, V4a-d, V6, V7, V8, V9, C107, C108. One BLOCKER under C100:
  `hardsoil_gammap_initial_index1.dat` (initial state via node_dof, not run by build_safe) lacked
  `global_element_dof_from_node_dof -yes` and failed its targets under the -no default (rc 1). Fixed by adding
  the record (output = baseline, gammap_initial 0.000776285656513) and pinned with run_ok in build_safe.
  Advisories (not fixed): mstrain_druckprag node_rhside <=2e-5 in near-zero slots; post_data missing index -> 0
  without a warning; plasti_apply2 disx/eptxx 4.66e-6 -> 0 (still = Pro references); the old text of
  mdiprisco_hist (without -yes) aborts at plasti.cc:485 (assert tmp!=0.) instead of a clean error.
- L-C110 (2026-11-05, owner, verbatim): "converjamos. Registramos la suma" -> C110 resolved inline (823d58b): signed sums in
  post_node_rhside_fixed/_free; ratio measured on Pro = max over principal dofs of mean|free|/mean|fixed| (mean free
  if fixed < 1e-10 or method -free); post_node_rhside_fixed_free() re-run after the post-solve pass (also resolves
  C63, slide4 = Pro). RED: old binary ratio 1.80/11.49 vs Pro 12.944/1.2388; GREEN 12.943/1.2388. Tests
  post_node_rhside_ratio1 (new), _method1 (updated). Suite 747 OK. Iteration counts change in the 2
  control_timestep_iterations_automatic inputs (default1 843->723, Pro 654; min_max_wished1 944->968, Pro 682).
  Pending note: the manual default dof_type -materi_velocity (6.951) not applied (resolved, L-T2-close).
- 2026-11-05 C110 independent verifier (fresh mktemp copies 155e5d9 / 823d58b): PASS. Suite 747 OK; 687-input
  corpus: rc identical, beyond the post_node_rhside records only the 2 control_timestep_iterations_automatic inputs
  differ (dbmeth20 = parallel-runner artefact, identical in isolation). Probes = Pro (dof_type filter honoured).
  Advisories documented: Pro accepts method -yes (GNU rejects); no re-run when the post-solve pass is skipped.
- L-T3 (2026-11-06) owner: "adelante con T3".
- L-T3-map (explorer, read-only): only hypo6.dat (niemunis_visco) and hypo11.dat
  (wolffersdorff_niemunis) of the 29 appear in Pro tests; GNU registers but never reads
  control_materi_plasti_hypo_substepping and group_materi_plasti_visco_exponential_*;
  GNU umat_ (umat.c) takes props from group_user_data; no PARDISO (so.cc dgesv).
- L-T3B: aliases in db_number; umat_parameters -> GROUP_USER_DATA (one GNU user-data
  record, documented). New global MATERI_PLASTI_MAXIMUM_ITERATIONS, group record wins.
  Pro probes (plasti_apply1 square): default -3.39353 = GNU; cap 1 -21.9773 vs GNU -21.978;
  Pro accepts 0 (GNU error, advisory). Manual default 20000 vs GNU 1000 -> BACKLOG C112.
  The old manual-user page described the Pro name with the group syntax; rewritten.
- L-C112 (2026-11-06) owner: "vamos a 20000. Es un límite que si hay convergencia más rápida no llega a cumplirse y evitarmos cortes antes de convergencia".
  MAX_PLASTI_ITER 20000 (membrane loop keeps MAX_ITER 1000). Impact (instrumented corpus, 931 inputs: test-2026, test-2014, Professional test/other, validation, incremental_driver): only sfnet `examp10` and `examp16` pass 1000 iterations, and both reach 20000 (points that do not converge); every other input is unchanged. `examp10` moves 3e-7 in the displacements; `examp16` node 26 sigyy -218129 -> -237309 (Professional -271692; summed |GNU - Professional| over node_dof 3.61e6 -> 3.44e6), run time 1 s -> 3 s. Test `materi_plasti_maximum_iterations_default1` (examp16 copy, node 26 target).
- L-T3C1: Pro takes `materi_plasti_hypo_substepping -no` without index (manual header shows one); the
  control record wins. GNU: thread-local switch skips the hypo.c accuracy halving (failure halving kept),
  minimum e in sigma_() and after the Masin kernel. Pro probes match to 2e-5 (hypo1, hypo7). Masin
  ignores the substepping switch in GNU, Pro moves hypo7 by 0.0036 -> DIVERGENCES HYPO-SUBSTEPPING-MASIN
  (class B).
- 2026-11-06 C1 independent verifier (Risk item 4, fresh copies ff42b7a / 6745219): PASS. TLS flags set before
  every read (hypo_ and sigma_ have one caller each; initial-exec TLS on both sides); 8x4 two-group probe
  with options_processors 1/4/8: group isolation exact, 1-vs-N differences 1.6e-8 = baseline threading
  round-off; 25 hypo inputs byte-identical to baseline without the records (time_calculation excluded);
  suite 757/0; rejected inputs rejected as Pro. Advisory: minimum_void_ratio1 hyhis0 tolerance 1e-5 pins
  GNU (Pro 1.2e-5 away).
- L-T3C2: mpc wall = node of the element in any active mpc_node_number (positive entries);
  plasti_on_boundary 2 selects group_materi_plasti_mpc_factor (plasti.cc, stress.cc cutoff); bounda wins
  (Pro: bounda+mpc = bounda alone). Pro probes (two quad4): locality, -no and factor-only match. Pro reduces
  phi AND c (explicit group phi,c x2/3 reproduces the Pro wall value -1.87260 exactly); GNU reduces phi only
  -> C114 (also bounda; direct full MC ignores the factor, plane cutoff default 1 vs 2/3).
- L-C114: manual and Pro agree (no owner decision, S2). Pro rule: on a wall the _wall record if given (no factor
  on top), else phi, phi_flow and c times the factor (default 2/3) - MC, hardening/softening, DP, direct full and
  plane. GNU fixed in plasti.cc (5 sites) and stress.cc direct_wall_factor(). Direct modes = Pro to all digits;
  incremental within 2.5e-3 relative (pre-existing MC gap). RED on a92f8bc build. Corpus wall inputs: only
  Pro other/mohr_coul_direct8 (identical to Pro) and the heat_exchanger_pile timeouts.
- 2026-11-06 C114 independent verifier (fresh a92f8bc / 397a3cb): PASS. Own probes re-derive the Pro rule (MC, DP,
  hardening, both direct modes, factors 0.5/0.6, _wall without factor on top); wall = explicit group bit-identical.
  75 plasticity inputs: only the 6 wall tests change, all toward Pro; no-wall inputs identical. Suite 765/0.
  Advisories: bounda detection broader than Pro + missing wall laws + single factor -> BACKLOG C115 / DIVERGENCES
  WALL-DETECTION; mpc page precision claim and bounda dev page corrected. Two advisories rejected with evidence:
  the C114 hash was already recorded (65771f0) and the manual numbering 6.687/6.688 is right.
- L-T3C3: phimob Pro probes: degrees (manual says both radians and degrees), root of the manual f with tension
  positive, root closest to 0 (5.3956 not 90), mohr_coul and mohr_coul_direct only. The dof label "phimob"
  collides with the GNU post_calcul operator PHIMOB (db_number finds the operator first): the label reuses
  -PHIMOB. residual_stiffness: Pro gives identical results and iterations with factor 1 / 0.1 / 0.01
  (plasti_apply1 square, examp16 180 iterations) -> owner decision (manual vs Pro).
- L-T3C3R (2026-11-06, user): "si, me parece bien" (to: implement residual_stiffness as the manual, matrix only,
  default 1, divergence with Pro). Implemented per integration point (plastic = yield test started plastic
  iterations or direct cut in the current update), C' = V + factor (C - V). Probes: 50 iterations disx
  0.0182467717 with and without factor 0.1; 1 iteration sigxx 0.0691216 (none, = Pro) vs -0.0101589 (0.1);
  Pro identical with the record and accepts 0 (GNU error).
- L-T3D1 (2026-11-06): elasti_c / c_direction. Pro probes: tensor strains (isotropic C = young / poisson, also
  in shear), additive with young / poisson, axes as rows of the rotation and normalized, membrane refused, zero
  axis degenerate (GNU error). Block placed after every elastic law (several branches clear C). GNU = Pro 1e-9.
- L-T3D2 (2026-11-06, user): "adelante con D2". Pro probes (1D visco bar): named model = visco_power, unused
  name elasto-plastic, visco_power wins over names, visco_apply -no gates them, full list = visco_power. GNU
  selects per integration point at the first yield test (visco_power_select). Pro reads
  visco_exponential_name as an INTEGER and rejects names ("Integer values are expected"), so the exponential
  twins contradict the manual 6.745 -> owner question. C116 added (tension + vonmises group differs from Pro).
- L-T3D2x (2026-11-07, user): "si, hazlo como en el manual y registra la divergencia con Pro". Exponential names
  share visco_named_select (visco_named 2); _limit (also registered, unread) replaces EPS_VISCO, Pro honours it
  (limit 1: -0.453047, GNU -0.453174). Named targets = GNU plain record (Pro refuses names). C117: plain
  exponential with the limit inactive differs from Pro (alpha 1, 0.01; not the iteration cap).
- L-T3D3 (2026-11-07, user): "adelante con D3". Pro other/hypo11 (wolffersdorff + niemunis 0.2 5.): Pro sigyy
  -0.976887 (target -0.976926), without the niemunis record -1.019954; GNU without it -1.15668 because the
  void_ratio_linear alias switched on the pressure dependent correction (C118). Fixed: GNU -1.0198992 (-yes),
  -1.0199336 (-no, Pro -1.019924). The manual gives b2 and f_d tilde only by reference to Niemunis 2003
  (eq. 4.211, 4.221), not reproducible from the manual -> owner question.
- L-T3D3b (2026-11-07, user): "te he dejado la tesis en external-downloads" (../external-downloads/
  niemunis-habil-2003.pdf). b2 eq. 4.183, f_d_bar eq. 4.221, interpolation 4.223. b2 in L only gave -64.5;
  with N updated through the unchanged flow rule (thesis below 4.183) GNU = Pro within 6e-5 for nu 0.2 / 0.365 /
  0.3822 / 0.49; dense states follow the Pro history to 0.1 % up to the collapse.
- L-T3D4 (2026-11-07, user): "adelante con D4". Manual 2.2.2 sentence on the Borja Tamagnini C_ijkl is cut
  ("is specified in" + nothing); only reference [1] Borja, Tamagnini 1998. Pro probes (hex8, homogeneous,
  G0 alpha khat pr; needs materi_strain_elasti): incremental and step dependent; isotropic tangent
  K = pr/khat max(1, exp(-eps_v/khat - 1)) at the start of the step (exact for dt 0.25 / 0.05 / 1, khat 0.01 /
  0.02, pr 5 / 10; extension stays linear pr/khat); pure shear G = G0 + alpha/khat (independent of pr);
  combined compression + shear couples (p -15.014062, G 311.25 in step 2), form not identified -> owner question.
- L-T3D4b (2026-11-07, user): "ya lo tienes en external-downloads" (../external-downloads/
  1-s2.0-S0045782597001412-main.pdf). Paper eq. 3.4-3.12: hyperelastic, P = P0 exp(Omega) (1 + 3 alpha eps_s^2 /
  (2 khat)), Q = 3 (mu0 - alpha P0 exp Omega) eps_s, Omega = -(eps_v - eps_v0)/khat; path independent, and
  at zero elastic strain P = P0. Pro is step dependent, zero stress at zero strain, shear G0 + alpha/khat
  (paper mu0 + alpha pr): Pro does not follow the paper -> owner question (paper vs Pro vs discard).
- L-T3D4c (2026-11-07, user): "1. El paper, y documentamos divergencia con Pro". Analytic branch in
  hyperela.cc (exact tangent), stress difference as the other hyperelastic laws; GNU = analytic to 1e-7.
- L-T3D5 (2026-11-07, user): "adelante con D5". Pro probes: components zeroed after plasticity, no
  condensation, axes normalized, zero axis runs; full tangent (5 / 30 iterations to converge on the manual 9.1
  layer). GNU projects stress and ddsdde after set_stress (materi.cc), automatic frame from ELEMENT_NORMAL;
  identical to Pro (converged). Pro needs global_normal -yes (GNU unknown, C14).
- L-T3E (2026-11-07, user): "vamos con la familia E", then "1, adelante" (to: low-cost hooks and pardiso now,
  diprisco_density now, niemunis_visco as a separate task at the end of T3). Pro: young_user -yes "Error: routine
  user_young is not programmed.", history_variable_user -yes "Error: routine user_materi_history_variable is not
  programmed." (Pro init record materi_history_variable); pardiso_decompose accepted.
- L-T3E2 (2026-11-07, user): "hacemos 1" (to: linear interpolation with the current void ratio, starting from
  history variable 12 and evolving with the volumetric strain, clamped to [ed, el]; class D). Pro probes: only
  the first set is used, el > ed required. GNU tests density1-4 (e 0.75 = averaged plain model exactly; e 0.9
  loose; e 0.5 dense; oedometer e analytic), error1 (Pro message), error2. Plain diprisco GNU vs Pro 0.2 % (C119).
- L-T3E3 (2026-11-07, user): "adelante" (E3), then "vamos con el recomendado" (to: b2 < 0 is an error per the
  manual, class D). Explorer: Jerman-Masin UMAT is not the Niemunis law; flow rule sign from habil eq. 4.73.
  Pro probes: hypo6 is pure relaxation; phi / nu / betaR do not affect isotropic relaxation; IGS required;
  OCR applied once in the first step of the index; b2 < 0 used. GNU within 5e-5 relative of Pro on 20 probes,
  Python reference equal to 1e-8. Tests visco1-4, ocr1-3, error1-2.
- 2026-11-07 end-of-T3 independent verifier (baseline 7101205, candidate f427646): PASS, no blockers. Suite
  833 OK; 224-input plasticity/hypo differential: only intended changes (C112 cap, C118 alias fix, T3
  keywords now accepted) and wall-clock lines; D2/D2x probes base = cand without T3 records. Advisories:
  unstable tension + vonmises visco (also at baseline, Pro NaN; C116), hypo10 own-target gap (new C121).
  S5 refresh: 131/320 of the 2026-11-01 list accepted, missing 189 = 63 GiD + 126 (KEYWORD-GAP-REPORT).
  T3 closed. Next: T4 W5 geometry, mesh tools, mpc, dependency.
- 2026-11-08 Owner: "Vamos con T4". L-T4-map (explorer): 69 names; aliases in db_number (database.cc ~9860);
  geometry entities in geometry.cc:geometry (new entities also in geometry.cc:920 lists, change.cc:81-90,
  check.cc:341); dependency in group.cc:get_group_data 458-555 (always "use", exact group match, no -all);
  refine_l.cc:126-139 hard-codes not=all / only=any; mpc options in mpc.cc:703-916 (MPC_EPS_ISO :76);
  control_mesh_truss_distribute_mpc and mpc_geometry are registered but never consumed; no geometry_moving,
  PRISM12/18 in GNU. Pro examples: validation_14_mesh (delete_direct), other/prism6, other/map1,
  other/safety2 (extrude_direction), other/extrude1 (element_group_new: layout differs from manual 6.184),
  other/truss10 (_air), other/geometry_moving1.
- 2026-11-08 U1: Pro accepts control_mesh_delete_geometry_element_group and geometry_triangle_eps_iso, rejects the
  manual plural mesh_interface_triangle_coordinates (DIVERGENCES INTERFACE-TRIANGLE-PLURAL, class C).
- 2026-11-08 U3 differential (baseline d6e9b27, 25 inputs with dependency_item or print_group_data): only distri3
  changes (print_group_data now shows the control_distribute values, as the Professional: ~100 instead of the raw 100000);
  examp24 and foundation_consolidation_settlement time out at 300 s in both binaries (rerun with 1800 s).
- 2026-11-08 U3: Pro probes (two-element square, mohr_coul radians): -all only for groups with the record, -multiply,
  number (row = number when one row per value), -tangent/-sinus/-cosinus, one row broadcast to every value, rows != 1
  and != nvalue error, geometry = all element nodes inside. All equal to GNU at the same evaluation time. C122:
  Pro evaluates -time_current at t_start + dt/2 (GNU end of step) and integrates total_linear per step (GNU E(t) eps,
  analytic); owner decision pending. print_group_data now shows get_group_data values (manual 6.406 sentence).
- 2026-11-08 U2: Pro probes (print_node_geometry_present, 2D 3x3 quad4 and 3D 2x2x2 hex8): all options/entities agree.
  Pro facts: quadrilateral eps = triangle rule; mpc -yes slave only, -no no restriction; element_geometry -any/-all/-only;
  cylinder angle refs z: x->y, x: y->z, y: x->z; self exclusion ignored by Pro; oblique cylinder without start vector
  selects nothing in Pro; Pro start-vector message copied. Pro defect on the plane through the first cylinder axis
  point (class D CYLINDER-PART-START-PLANE, S2 rule, owner confirmation pending).
- 2026-11-08 U4 (route delegated, context backstop): one shared element_in_geometry() (geometry.cc, -all/-any/-average)
  used by element_geometry_present, dependency_geometry, delete_geom and refine_locally. Pro probes (4x4 quad4 with a
  quadrilateral over element 1, circle 2 2 r 1.45; two quad4 with a dependency point; eight tria3 refined with -nothing):
  geometry_method honoured by element_geometry_present (-any 1 2 5 6, -average 1) and dependency_geometry (-any 300/300,
  -average 300/2000), ignored by control_mesh_delete_geometry, refine not/only and area_element_group; delete -any keeps
  3 4 7 8 9-16, -average deletes 1; _direct -yes removes completely in the first step (factor ignored); refine only default
  -all = manual (16/14; GNU was -any 32/25), only -any 32/25, not default 28/22, not -any 8/9, minimal_size 1.5 8/9,
  1.4142 32/25: GNU = Pro on all. Pro circle delete tests the edge only and ignores _projection_type and
  geometry_projection_type: GNU keeps the 6.166 filled default, -project_exact per 6.173 (class D DELETE-PROJECTION-TYPE,
  owner decision C123). Invalid values: GNU errors, Pro ignores / asserts (class C U4-METHOD-VALUES). C124: GNU does not
  print element_delete_factor (Pro does). Word -project_exact registered (was rejected). 24 tests, suite 917 OK / 0 FALLO,
  audit REVIEW 42, 26 warnings = baseline. Differential vs d467efa (29 inputs using delete / refine_locally /
  print_element_geometry_present / dependency_geometry, suites + Pro test dir): every .dbs identical (time_calculation
  ignored), same rc; tutorial_3 and slope_nonlocal_refine hit the 60 s timeout in both, validation_14_mesh fails to parse
  in both (change_dataitem_time_until_data), so its _direct example is not reachable.
- 2026-11-08 U3 long differential (1800 s): examp24 identical .dbs; foundation_consolidation_settlement times out in
  both binaries (one diagram row, one value: unchanged code path by inspection); left to the end-of-T4 verifier (omen).
- 2026-11-08 owner on C122 (verbatim): "yo creo que mantenemos GNU, es fácilmente documentable y reproducible, y no
  nos ha dado fallos. De todas formas, documéntalos". Recorded: DIVERGENCES DEPENDENCY-TIME moved B -> D (B 55, D 21),
  SEGUIMIENTO Descarte row, user/developer dependency_item Differences, BACKLOG C122 decided (59 open). Next: C123.
- 2026-11-08 owner on C123 (verbatim): "supongo que lo correcto es el círculo relleno. Sería muy raro y bizarro eliminar
  un borde. Generaría superficies no conexas". Recorded: DIVERGENCES DELETE-PROJECTION-TYPE decided (class D), SEGUIMIENTO
  Descarte, element_geometry_selection user/developer Differences, BACKLOG C123 decided (58 open). Code unchanged.
- 2026-11-08 owner confirms CYLINDER-PART-START-PLANE class D (verbatim): "si". Recorded in DIVERGENCES, SEGUIMIENTO
  Descarte, geometry_options_2026 user/developer Differences. Next: T4 U5 (mesh operations).
- 2026-11-08 L-T4-U5-map (explorer): none of the 15 names registered. Pro examples: other/extrude1
  (_element_group_new, matches manual 6.184: old groups, then n_old x n_layers new groups), other/safety2
  (_extrude_direction -y), other/prism6 (_split_element_to -prism6); the rest have none. GNU mesh_extrude
  (mesh.cc:852) keeps the 2D nodes at z=0 as an extra bottom layer and ignores _n for tria3/quad4/bar2. Probe
  (build of 87436f7+): Pro extrude2 3 el / 12 nodes, interface_bar2_hex8 6 el / 24, extrude1 45 el / 119; GNU
  interface_bar2_hex8 7 el (C50 symptom), extrude2 stops "prism6 not available for force_element_edge",
  extrude1 rejects _element_group_new. mesh_remove (delete.cc:894) lacks the -method3 branch its comment names.
- 2026-11-08 U5a (route inline): mesh_extrude(icontrol) rewritten to manual 6.179 / 6.185 with the Professional
  numbering (plane p of node i = p*n+i, layer l of element e = l*m+e, from db_highest_index), node orders (quad4 ->
  hex8 a b a' b' c d c' d', so interface quad4 sides stay faces) and unused-node deletion; extrude() runs once at its
  own control index (also without control_timestep: Pro other/extrude1). -y: the 2D y becomes z (Pro probe). input.cc
  reads negative numbers of _element_group_new. One coordinate kept as a GNU extension (documented GNU usage; Pro
  ignores it). Pro = GNU .dbs mesh: extrude1, interface_bar2_hex8 (+ node_dof 7.6e-12, C50 mesh resolved), probes
  layers / direction / group_new / spring1-2. Differential vs HEAD d08a5d4 (8 inputs): results equal by coordinates,
  numbering compact (extrude_quad9_shared_nodes1 targets 83/66/71 -> Pro 45/30/33); extrude2 now aborts in point_el
  (no prism6, pre-existing, C125, to fix in U5b) instead of the force_element_edge error. RED: base binary rejects
  the new records / leaves extrude_layers1 unextruded (no timestep). Risk: item 3 (control_mesh_extrude existing
  inputs change numbering and layer convention; covered by the end-of-T4 verifier).
- 2026-11-08 U5b (route inline): split.cc _element_from filter, _element_to -prism6 (2 prisms, Pro order), tria6 ->
  4 tria3 (manual 6.240, missing before); mesh_convert_tria6_tria3 (corners n0 n2 n5, unused nodes deleted, as Pro);
  point_el prism6 via the prism15 Newton (C125 part). Pro = GNU: other/prism6 (targets), probes split_from2 /
  split_mixed / conv1. Pro ignores an invalid _to (GNU error, class C SPLIT-INPUT-CHECKS). Differential vs bccfbba,
  33 inputs: only prism6 changes (now rc 0). extrude2 back to the force_element_edge prism6 error (C125 open part).
  RED: base rejects the records / leaves tria6 unsplit. Suite 941/0, audit REVIEW 42. Risk: item 3 (mesh_split and
  point_el shared; no existing output changed in the differential).
- 2026-11-08 U5c (route inline): offsets in the top.cc control_mesh_renumber block (explicit records: a record number
  in a variable raised the audit to REVIEW 44), mesh_duplicate_element_group (delete.cc), control_mesh_remove_geometry
  through element_in_geometry -all. Pro probes (three quad4): same elements / groups / element_geometry. Differential vs
  b88adc5, 6 inputs: only the new tests change (base rc 1 = RED). Suite 953/0 before the audit rewrite. C126 found.
  Risk: none beyond item 3 of mesh_remove (shared; remove_method1_1 pins the old -method1 result).
- 2026-11-08 U5d (route inline): generate_beam_truss separate -yes = every later geometry node instead of NODE_NODE
  (Pro: all pairs in node order); mesh_correct_interfaces (mesh.cc) with the Pro permutations (quad4 xi -> n1 n3 n2 n4,
  hex8 xi -> n2 n4 n6 n8 n1 n3 n5 n7, hex8 eta -> eta/zeta swap, ties kept). 14 Pro probes equal. Differential vs
  394ac2a (94 inputs with generate truss/beam or group_interface): identical except the new tests (base rc 1 = RED);
  foundation_consolidation_settlement times out in both. Suite 971/0, audit REVIEW 42. Risk: item 3 (generate_beam_truss
  shared; generate_truss_line1 pins the old path).
- 2026-11-08 owner on U7 (verbatim): "implementamos la feature completa". U7 = implement control_mesh_truss_distribute_mpc
  (6.245: truss nodes tied with mpc to the isoparametric elements they run through; _exact redistribution, its minimal
  length / _connect, _dof, _element_group_truss / _isoparametric, _geometry_isoparametric) and mpc_geometry_method /
  _tolerance. Pro examples: other/truss10-15, mpc2, earthquake_1-3. Starts after U6 (shared mpc.cc).
- 2026-11-08 U6 (delegated writer): mpc.cc mpc_element_group_generate split into mpc_element_group_index (per index)
  with _always -no (skip nodes of element_group_1 elements), _closest (Euclidean closest point for bar/quad/hex,
  Professional clamped barycentrics for tria/tet), _coord_geometry (last coordinate >= lowest element_group_1 node),
  _eps_iso (point_el tolerance, default 1.e-4 kept), _keep (records outside the regenerated range, internal record
  mpc_element_group_kept) and control_mpc_element_group -yes (fingerprint xor TIME_CURRENT hash). No Professional
  example uses the records; ~60 Pro 02-08-2026 probes: equal for always -yes, closest on rectangles / hex8 / tria3,
  coord_geometry with the lowest node on the geometry (2D y, 3D z), eps_iso 1.e-3 / 1.e-4 tie sets, control -yes.
  Pro contradicts the manual (GNU follows the manual, S2; class D rows MPC-ELEMENT-GROUP-ALWAYS / -COORD-GEOMETRY /
  -EPS-ISO / -KEEP, CONTROL-MPC-ELEMENT-GROUP): ignores always -no and keep, coord_geometry only over element_group_1
  nodes on the geometry (no effect without it), default eps ~1.e-3 and eps >= 1.e-2 ties nothing outside, re-creates
  every timestep with control -no. Class C MPC-ELEMENT-GROUP-CLOSEST (distorted quad4 approximation, tet4 ignores
  the fourth corner); C127 found (extrapolated factors / point_el corner misses for outside nodes within eps).
  RED: the base binary rejects the 12 new-record tests (rc 1); base1 and control_mpc_element_group_no1 pin the old
  behaviour. Differential vs the U5d binary, 35 inputs with mpc records: 29 SAME, 5 rc 1 in both, truss13 executable
  path only. Suite 999/0, runs 226/226, audit REVIEW 42. Not edited (outside the surface): the older pages
  manual-user/mpc_element_group.md:63 and manual-developer/mpc_element_group.md:33,65-70 still call the member skip
  "always -no" and list the options as pending. Risk: item 3 (mpc_element_group_generate shared by every
  mpc_element_group input; base1 and the differential pin the old path).
- 2026-11-08 U7 (delegated writer, S1-S5): mpc.cc mesh_truss_distribute_mpc() (called once from step_start after extrude at
  its own control index): -truss elements recreated after the highest element (old one deleted, as Pro), _exact split at the
  entries into a next isoparametric element (Cyrus-Beck clipping td_clip on the corner sides; bar/tria/quad/tet/hex/prism6,
  no interface groups), minimal length (default 1.e-11 x truss length, Pro probes) and _connect (default -yes, Pro), _air
  (centre rule), _dof, truss/iso group or geometry selectors paired by position (Pro error texts for both-given / unequal
  counts / two records), ties via point_el + mpc_weights_polish through the shared writer mpc_tie_put() (also used now by
  mpc_element_group_index), marker mpc_from_truss_distribute_mpc; air end of an exact piece tied at the exit point (Pro).
  mpc_geometry_generate() (from mpc_node_apply, fingerprint): -method0 (switch + tolerance), -method1, -method2, dof order
  as listed, marker mpc_from_geometry, no factor record (Pro); _dof and _switch (method0) required as the Pro. New word
  -method0. Verification: Pro examples truss10/11/12/13/15 and mpc2 copied to test-2026 (targets met; truss13 with the GNU
  element number 6 and 50 iterations, C128) and ~50 Pro probes: every test .dbs equals the Pro one after matching nodes by
  coordinates except the manual cases. Class D (GNU = manual, no owner decision recorded): TRUSS-DISTRIBUTE-AIR-DEFAULT
  (manual default -yes, Pro behaves as -no), -AIR-NO (Pro deletes a truss with any node outside; manual: centre),
  -AIR-EXIT (Pro -yes also splits at the exit into the air), -ISO-SELECTION (Pro ignores _isoparametric restrictions).
  Class C: TRUSS-DISTRIBUTE-NUMBERING (compact numbers; Pro leaves gaps), MPC-GEOMETRY-METHOD2-GROUP. _dof indexed as the
  Pro although the manual heading shows no index. Not implemented: -trussbeam generation (no record; Descarte). C128 (mpc
  tie convergence, Pro truss13 3.87/6.13 with 2 iterations, also with hand-written records before U7) and C129
  (bounda_time limited to 190 values: Pro earthquake_1-3 stop at input; they are the mpc_geometry examples) found.
  RED: the base 4bfbe19 binary fails 23 of 24 td_check cases (truss11 pins the old behaviour) and the 5 error tests.
  Differential vs 4bfbe19 (101 inputs with truss/mpc records): 87 SAME, 5 rc 1 in both, intended changes truss10/12/13/15
  and mpc2, executable path / time_at_start only otherwise. Suite 1052/0, runs 226/226, audit REVIEW 42. Risk: item 3
  (mpc_element_group_index now writes through mpc_tie_put; mpc_node_apply also generates mpc_geometry; the
  _element_group_truss record became a list; the suite and the differential pin the existing paths).
- 2026-11-08 U8 (delegated writer, S1-S4): new module geometry_moving.cc; control_geometry_moving -initialise (once at its
  control index, step_start) samples nspace^d points per element (manual 6.538; isoparametric cell midpoints, Duffy for
  simplices, nodes for other elements) against the entity with geometry() at ntime+1 moments (points moved back over the
  translation, scratch geometry record); records element_geometry_moving / _part / _in_time / _out_time (Pro names);
  geometry_moving_apply every step_start: delete factor ramps 1 -> 1-part between in_time[0] and [1] (Pro measured),
  0 deletes the element and its unused nodes; ELEMENT_RHSIDE_DELETE allocated before the parallel loop. Pro probes
  (other/geometry_moving1 with processors_maximum -no, ~40 runs): time conventions and ramp reproduced; Pro tests the NODES
  and ignores nspace -> class D GEOMETRY-MOVING-SAMPLING (manual followed), GEOMETRY-MOVING-PARAMETERS (Pro quad lengths),
  GEOMETRY-MOVING-THREADS (Pro abort), class C GEOMETRY-MOVING-CONTROL-VALUE; no owner decision recorded. Tests:
  geometry_moving_band1/_band2/_simplex1/_hex1/_ramp1/_ramp2/_error1-5 (analytic), geometry_moving1 Pro copy (not in the suite:
  35 min, GNU sigyy -0.3489 vs Pro target -0.3268 +- 0.02, rc 1; bicg breakdown, BACKLOG C130). RED: base rejects
  geometry_moving in all 11 tests. Differential vs HEAD (45 inputs with change_geometry / change_dataitem /
  control_mesh_delete_geometry, validation-suite + Pro test/): 35 SAME, 3 SAME-NODBS (rc 1 both), 4 DIFF-OUT executable path
  only, 3 rc 124 in both (timeout 300: slope_nonlocal_refine, foundation_consolidation_settlement, dam_building; output
  length only). Suite 1068/0, runs 226/226, audit REVIEW 42. Risk: item 3 (step_start calls the new hooks for every input;
  they return at once without the records, the differential pins it).
- 2026-11-09 owner confirms the four U7 class D rows (TRUSS-DISTRIBUTE-AIR-DEFAULT / -AIR-NO / -AIR-EXIT /
  -ISO-SELECTION: GNU follows the manual), verbatim: "si". Recorded in DIVERGENCES and the SEGUIMIENTO Descarte row.
- 2026-11-09 owner (verbatim "si"): U8 class D confirmed (GEOMETRY-MOVING-SAMPLING / -PARAMETERS) and geometry_moving1
  retargeted to the GNU value: full run (35 min, build of 31c5432) sigyy -0.348912581336, target tolerance 1.e-6 (the
  Professional -0.326823 kept in the header). Still out of build_safe (C130 run time).
- 2026-11-09 owner, U9: implement the new elements instead of discarding control_mesh_extrude_element and
  control_mesh_convert_quad9_quad6 (prism12, prism18, continuum quad6), verbatim: "A mi me parece que es interesante
  añadir los elementos. No es algo especialmente complejo, ya que tenemos toda la infraestructura y pueden derivarse de
  los que ya existen, y dota de escalabilidad al proyecto a futuro". The derivation claim is being verified by a
  read-only map (route: delegated explorer, mapping trigger).
- 2026-11-10 map result: all three shape families are products of existing functions; continuum quad6 can share the
  QUAD6 enum because interface status comes from GROUP_INTERFACE (elem.cc:515); MNOL 27 is enough. Manual gaps go to
  Professional probes: integration rules and point order, quad6 node order after conversion. First prism15 commit
  443fdcd shipped a wrong triangle table without tests, so each element lands with an analytic integration test.
- 2026-11-10 owner, VTK scope, verbatim: "si, la incluimos. Si vtk no lo soporta directamente convertimos en la
  escritura del vtk a un tipo soportado". U9 split into U9a, U9b, U9c.
- 2026-11-10 U9a parent review: the writer wrote quad6 as VTK cell 31 (VTK_QUADRATIC_LINEAR_WEDGE) and pinned the same
  value in its test. Fixed to 30 (VTK_QUADRATIC_LINEAR_QUAD). Checked with VTK 9.7.2 reading the files: quad6 and prism15
  cells have the right class, node count and mid-edge positions. Differential of 45 inputs against 408d8a1 identical.
- 2026-11-10 owner, QUAD6-3D (Professional runs quad6 in 3D, the manual lists it for 2D only, GNU refuses it),
  verbatim: "si, confirmado". Class D stays in DIVERGENCES-PROFESSIONAL.md.
- 2026-11-10 U9b writer result: prism18 = tria6_shape x 3-node line, layered tria6 node order, 21 points for the
  default and -maximal (= Professional to 9e-14). Professional -minimal (8 points) gives an indefinite stiffness
  (6 negative eigenvalues, pardiso -4); GNU keeps 21 points as C82 did for prism15, owner question C132. Recorded
  divergences: PRISM18-POST-POINT D (Professional interpolates linearly), EXTRUDE-ELEMENT-MIXED D (Professional
  gives prism12 for tria6+quad4), EXTRUDE-ELEMENT-NAME C (unknown words are an error), EXTRUDE-ELEMENT-PRISM12 B
  (refused until U9c), VTK-PRISM18 A (type 32; Professional writes a linear wedge). Parent review: VTK 9.7.2 reads
  the cell as vtkBiQuadraticQuadraticWedge with every mid-edge and face node in place; differential of 278 inputs
  against 76b0759 identical; stale lines in prism15.md and control_mesh_extrude.md (user and developer) fixed by
  the parent. Follow-ups: C133 (prism18 face loads, refine, split, mpc polish, geometry_moving), print_vt.cc
  length loop still lacks -TRIA6 and -TET4 (pre-existing).
- 2026-11-10 owner, PRISM18-MINIMAL (keep 21 points for prism18 -minimal, as C82 for prism15), verbatim: "si".
  Moved to class A in DIVERGENCES-PROFESSIONAL.md; BACKLOG C132 closed.
- 2026-11-10 U9c writer result: prism12 = tria6_shape x 2-node line, tria6 bottom and top layers, 14 points
  (7-point tria6 x 2 Gauss) for every group_integration_points value. Professional runs only -minimal (8 points,
  indefinite: 3 negative eigenvalues even undistorted, pardiso -4); the other rules stop with
  element_intpnt_iso_coord. Rows: PRISM12-MINIMAL A (C82/C132 rule), PRISM12-DEFAULT D, PRISM12-POST-POINT D
  (Professional interpolates linearly), PRISM12-INTERFACE B (refused, manual 6.624/6.158 list it),
  EXTRUDE-ELEMENT-PRISM12-QUAD9 D (Professional writes prism18 next to quad9; GNU follows 6.181 and the mesh is not
  conforming on hex27 faces), VTK-PRISM12 A (type 31; Professional drops the cell). EXTRUDE-ELEMENT-PRISM12 B
  resolved. Parent review: VTK 9.7.2 reads vtkQuadraticLinearWedge, GNU map equal to 2e-16; differential 289 of 290
  identical (the one difference is the deleted U9b refusal test); no attribution lines.
- 2026-11-11 owner, EXTRUDE-ELEMENT-PRISM12-QUAD9, verbatim: "la malla debe ser conforme. Creo que lo suyo es generar
  un prism 18, pero me plantea un problema, porque al otro lado del prism18 que hemos generado se estará generando un
  prism12 que volverá a no casar con el prism18 que acabamos de hacer. ¿Se podría generar un elemento de interfaz que
  casase las dos caras?". Probe: the Professional turns every tria6 into prism18 when the mesh has a quad9 (also a
  tria6 sharing only a corner node), so no prism12/prism18 face occurs and no transition element is needed. Proposal
  "-prism12 pasa a prism18 en toda la malla si hay algún quad9, con aviso" accepted, verbatim: "si". Implemented
  inline (RED: prism12 and no warning; GREEN, suite 1119 OK, 0 FALLO); the class D row was removed.
- 2026-11-11 owner, on PRISM12-DEFAULT (GNU 14 points of manual 6.623; Professional stops) and PRISM12-POST-POINT
  (GNU analytic; Professional corner-only interpolation), verbatim: "confirmamos clase D". Both rows stay class D.
- 2026-11-11 owner, on control_mesh_map (implement vs class B), verbatim: "si acelera el cálculo creo que merece la pena
  implementarlo. Es trabajo mediano". New task U9d. Existing map.cc is the remesh data map (map_element / map_node),
  a candidate helper for the linear -> quadratic interpolation.
- 2026-11-11 U9d writer result: new module mesh_map.cc (step_start after extrude, once per control index). Professional
  probes (.dbs of every stage): -hex27 -> -hex8 of the corners, mid nodes kept unused, element_map kept after the map
  back, each map toggles, mid node dofs after the map back = trilinear interpolation of all corner dofs (4.5e-12),
  -hex20 -> -hex8 1-8, prism18 / prism12 / quad9 / tria6 unchanged, a map at a timestep index runs before it,
  global_element_dof_apply not checked; the Professional maps -tet10 to 1 4 6 10 and -prism15 to 1 3 3 4 5 6. The GNU
  deletes unused nodes on every mesh change (nod_nod, node.cc outside the surface), so the mid node coordinates are
  kept in the GNU record element_map_node and the nodes created again on the map back. create_element / create_node
  are not usable in place (DATA_ITEM_SIZE copy buffer overflows on element_matrix_* of a solved element; db_dbl
  pointers invalidated by a PUT beyond the highest index), so mesh_map has its own in-place element / node PUTs;
  a shorter ELEMENT PUT needs a delete first (db_len bisection on the LONG_MIN end mark). GNU = Professional on a hex27
  column (1.3e-12 linear, 2.5e-12 after the map back, 3.7e-9 final). The map result is not the plain quadratic
  result (15 % on a mid node of the coarse column, = Professional): C134 owner question. Pre-existing findings:
  other/map1 hex27 ring -disz 1245 GNU -1.7795e-05 vs -1.9017e-05 also without the map (C135); -no with hex27 differs
  (C136); GNU control_print_database / control_print_gid ignore or refuse -separate_index (pre-existing, tests stop
  at the checked stage). Rows: MESH-MAP-MID-NODES B, MESH-MAP-HEX20 B, MESH-MAP-TET10 D, MESH-MAP-PRISM15 D,
  C135-MAP1-RING B, C136-GEDA-NO-HEX27 B. RED: base a565fc5 "I do not know : control_mesh_map"; GREEN: 14 tests,
  29 checks, suite 1148 OK, 0 FALLO (226/226 runs); differential of 144 mesh inputs against a565fc5 identical except
  time_calculation. Commit b3f4746. README indexes (S4) not edited (outside the writer surface).
- 2026-11-11 owner, on C134 (quadratic step keeps the interpolated linear stresses, = Professional), verbatim: "la contra
  que dices es relativa, porque tampoco tenemos datos para hacer una interpolación cuadrática. Reiniciar el cálculo en
  cuadrático también lo puede hacer el usuario "a mano". Creo que es mejor seguir al Pro en este caso". C134 closed.
- 2026-11-11 owner, on MESH-MAP-TET10 / MESH-MAP-PRISM15 (GNU true corners; Professional wrong corners), verbatim:
  "si, clase D para los dos".
- 2026-11-11 owner, on MESH-MAP-MID-NODES (keep element_map_node or change nod_nod), verbatim: "me parece bien
  mantenerlo así". U9d closed. Next: owner decision on control_mesh_remove_really (+_activate_all/_factor, _keep_geometry).
- 2026-11-11 owner, on the control_mesh_remove_really family (implement vs class B), verbatim: "me gustaría
  implementarla. Realmente no es muy costoso (creo) y: 1. mejora la convergencia; 2. mejora el flujo de trabajo".
  New task U9e.
- 2026-11-11 U9e writer result: new module mesh_remove.cc (mesh_remove moved from delete.cc, evaluated on every
  step of an index with control_mesh_remove). Professional probes (~60, 2D quad4 column, removed element on a
  coincident container of another group, container deleted, reactivation): default _really -yes; the Professional keeps
  the inactive element (element_remove_active -yes, element_empty -empty, cleaned) and its free nodes fully bounded;
  a step without control_mesh_remove changes nothing; factor / _activate_all of the reactivating index; mean over the
  element's active nodes (node_active_last) of all dofs, stresses and pore pressure times the factor, temperature,
  velocity, displacement unscaled; elements in order, a node set once, zero without active nodes; _keep_geometry of the
  current index (all nodes) followed by the index's _really; a reactivated element is not "removed once" any more;
  factor outside [0, 1]: "Error, CONTROL_MESH_REMOVE_REALLY_ACTIVATE_FACTOR factor should be between 0 and 1.".
  With global_element_dof_apply -yes the Professional starts the reactivated element from zero integration point
  stresses, so the factor does not change the result (identical for 0 / 0.5 / 1); the GNU mirrors it
  (ELEMENT_DOF_INITIALISED 1 over a zero ELEMENT_INTPNT_DOF) and equals the Professional on every node dof at every
  stage (5e-11). GNU representation as U9d: inactive element deleted and kept in element_remove_inactive / _group /
  _node, unused nodes deleted by nod_nod and created again (REMOVE-INACTIVE-RECORDS, B). Manual followed against the
  Professional (class D, owner confirmation pending): REMOVE-NO-TIMESTEP (Pro stops on a removal at an index without
  timestep, leaves a reactivated element inactive there), REMOVE-ACTIVATE-ALL (Pro also resets velocity and
  displacement of the active nodes), REMOVE-PRESSURE-GRAVITY (Pro factor * (mean - rho g.x_middle) + rho g.x_node);
  REMOVE-REALLY-WORD (C, GNU input error, Pro takes another word as -no). Pre-existing: -method1 matches node numbers
  (C126), so a geometric container with own nodes removes nothing in the GNU; global_element_dof_apply -no also
  differs on a quad4 column (C136 widened). RED: base 896978f "I don't know what to do with :
  control_mesh_remove_really", 32 of 34 new checks fail (control_mesh_remove_default1 pins the old deletion). GREEN:
  14 tests, 34 checks, suite 1182 OK, 0 FALLO (226/226 runs); differential of 93 inputs (53 remove / delete /
  keep / excavation / geometry_moving / element_empty / mesh_activate / map / renumber inputs + 40 random test-2026)
  against 896978f identical (geometry_moving1 times out in both). Commit 5f08f81. README indexes (S4) not edited
  (outside the writer surface).
- 2026-11-11 owner, on the U9e factor question (Professional starts a reactivated element from zero intpnt stresses,
  so _activate_factor has no effect), verbatim: "Lo que hace el profesional tiene sentido, pero se consigue poniendo
  un factor 0, así que si arrancamos desde los valores nodales cubrimos tanto el pro como el manual, ¿no?". Parent:
  only with factor 0 + _activate_all -yes (shared nodes keep their stresses), and the default factor 1 then differs
  from the Professional (class D, recipe documented). Owner, verbatim: "si". U9e reopened: start the reactivated
  element from the nodal values.
- 2026-11-11 U9e nodal start writer result (9c08983): the ELEMENT_DOF_INITIALISED 1 override is removed, so
  create_element's 0 makes the reactivated element interpolate all integration point dofs (stresses included) from the
  nodal values set on reactivation. Tests control_mesh_remove_start1 / _start05 (all nodes bounded, nodes 4 5 moved
  only in stage 10 by bounda_time, so no increment at 20: element_intpnt_dof stresses = element_intpnt_h x nodal start;
  factor 1 analytic -11.111 continued at nodes 4 7), _recipe_really1 / _recipe_chain1 (factor 0 + _activate_all -yes =
  the Professional default runs t1 / tc1, all 10 node dofs, 5e-11 and 1e-11). really1 / chain1 (factor 0.5) re-pinned
  to the GNU values (up to 1.07 / 1.31 from the Professional). Professional: its own factor 0 + _activate_all -yes run
  differs from its default run only by node 4 disx 8.06e-4 (REMOVE-ACTIVATE-ALL displacement reset); the GNU recipe
  reproduces the Professional default exactly. Recipe measured with an elastic material: in the GNU strains and
  plasticity history also start from the nodes. RED: the ac4bc9b binary fails start1, start05 and the start1 pin
  (3 of 43); GREEN 43 of 43, suite 1191 OK, 0 FALLO (226/226 runs); differential of the 93 U9e inputs against
  ac4bc9b identical (geometry_moving1 times out in both, common output prefix identical). Class D row
  REMOVE-REACTIVATION-START with the recipe, user / developer pages and SEGUIMIENTO updated.
- 2026-11-12 owner, on REMOVE-INACTIVE-RECORDS (GNU keeps inactive elements in its own records, = MESH-MAP-MID-NODES
  choice), verbatim: "si, clase B". Open: class D of REMOVE-NO-TIMESTEP / -ACTIVATE-ALL / -PRESSURE-GRAVITY; C126.
- 2026-11-12 owner, on REMOVE-NO-TIMESTEP / REMOVE-ACTIVATE-ALL / REMOVE-PRESSURE-GRAVITY class D, verbatim: "si".
  Open: C126 (-method1 geometric containment).
- 2026-11-13 owner, on C126 (implement -method1 geometric containment with -method3), verbatim: "si, me parece
  correcto". New task U9f.
- 2026-11-13 U9f writer result (d13d4e2): -method1 by position in mesh_remove.cc (reference elements collected per
  pass with a bounding-box prefilter, each candidate node located with point_el on the full node list of a reference
  element) and -method3 (node flags from the constrained node of the mpc_node_number records; METHOD3 word added).
  Professional probes (~60): union of references node by node (an element over two references and over a gap between
  them is removed: only nodes are tested), boundary inside, tolerance windows quad4 / hex8 [2e-4 in, 2.05e-4 out] in
  [-1, 1] coordinates (scales with the element), tria3 [0.825e-4 in, 0.9e-4 out] barycentric; GNU eps 2.02e-4 /
  0.85e-4 (tet / prism not measured, simplex value); element_group_0 elements are no references (-method1 2 2 removes
  nothing; the old node match removed one of two coincident elements); every node of a quadratic element tested; a
  curved quad9 reference is its straight corner element in the Professional, curved in the GNU (REMOVE-METHOD1-CURVED,
  D); -method3: masters do not count, groups stop with "Error detected for data item : control_mesh_remove";
  -method3 with _really -no / control_data_delete of an mpc / _keep_geometry = Professional states; the GNU runs with
  mpc records on nodes it deleted while inactive. Reactivation with a geometric container (p7) and with -method3 (w1):
  GNU recipe (factor 0 + _activate_all -yes) = Professional default to 5e-11 / 4e-11. RED: base 31c00c6 fails 22 of
  the 40 new checks (node match, -method3 unknown); GREEN 83 of 83 U9e + U9f checks, suite 1231 OK, 0 FALLO (226/226
  runs). Differential of 111 inputs (the 93 U9e inputs + every tracked input with control_mesh_remove) against
  31c00c6: identical except the tochnog_version build-date line (geometry_moving1 times out in both, common prefix
  identical); no existing input used non-node-sharing -method1 meshes. C126 closed in BACKLOG.
- 2026-11-13 owner, on the -method1 nodes-only test (Professional removes an element spanning a gap between containers),
  verbatim: "creo que es más riguroso seguir al manual, aunqeu más complicado. ¿crees que hay muchos casos en los que
  vayamos a notar la divergencia, o la limitación del Profesional es asumible?". Parent: visible only when an element
  is larger than a gap, a reentrant corner or a curvature of the reference region; proposed also testing the
  integration points of the element with the same locator (class D). Owner, verbatim: "si, añadámosla y lo
  documentamos". U9f reopened.
- 2026-11-13 U9f integration points writer result (728e997): -method1 also locates the integration points of the
  element (pol() h x node coordinates, as element_intpnt_set; an inactive element gets its ELEMENT record only for the
  pol() call; no integration rule or interface group: nodes only). Tests RED first: control_mesh_remove_method1_gap2
  (gap 1.5 < x < 1.7 under the integration point x = 1.577) and _lcorner1 (L of [0,2]x[0,0.5] and [0,0.5]x[0.5,2],
  convex quad4 with the integration point (1.577, 0.644) outside the reentrant corner) fail on d13d4e2 (2 of 87),
  GREEN 87 of 87; the Professional removes both (nodes only). Limits documented: _gap1 (gap 1 < x < 1.5 between the
  points) still removed; a quad9 bulging out of a convex reference element between its nodes is not constructible with
  an integration point outside (random search of 200000 valid quad9 layouts in a square found none), so no bulge test;
  a first L geometry (quad4 vertices 0.1/1.9) put no integration point outside either. Suite 1235 OK, 0 FALLO (226/226);
  differential of the 111 U9f inputs against d13d4e2 identical except the build-date line (geometry_moving1 times out
  in both). Class D row REMOVE-METHOD1-INTPNT, Descarte row, user / developer pages, C126 note.
- 2026-11-13 owner, on REMOVE-METHOD1-CURVED class D, verbatim: "confirmamos clase D". Tolerances kept as
  implemented (parent proposal, no objection). U9f closed; every U9 owner decision taken.
- 2026-11-13 end-of-T4 independent verifier (HEAD 3b2dd42 vs baseline cdcf986): S1, S2, S3, S5, S6 PASS; S4 one
  blocker: control_mesh_map and control_mesh_remove_really missing from both README indexes (fixed in the next
  commit, with the stale mpc_geometry "Pending" note of manual-developer/mpc_element_group.md). Suite 1235 OK / 0 FALLO,
  226/226. Invariant: 772 baseline inputs (test-2026 528, test-2014 244) run on both binaries, 758 identical, 14
  explained (time_at_start, gz, bicg_csr_threads1 thread noise, C91 enum shift, expected new-name lists, intended U5a
  extrude renumbering). Rejected inputs stop with the documented messages. Advisory: keyword probe accepts 68 new W5
  names (tracker says 69). Pending for T4 closure: KEYWORD-GAP-REPORT W5 refresh (S5).
- 2026-11-13 S5 refresh at the end of W5: probe on the build of 16f3721, 199 of the 320 names of the 2026-11-01
  list accepted (68 more than after W2/W3), missing 121 = 63 GiD + 58 pending (T5, T6 safety_*, 2 initialisation
  literals); none of the 58 is a T4 name. The tracker 69 vs observed 68 is not explained. KEYWORD-GAP-REPORT
  updated. T4 closed. Next: T5 W6 other physics.
- 2026-11-13 L-T5-map (explorer): 24 names + 2 initialisation literals (materi_history_variable,
  materi_strain_energy, input.cc:195-198) = 26 (task line says 25). No T5 name in the Professional test/ tree.
  Aliases: group_contact_spring_plasti_cohesion / _friction_automatic (conspr.cc:76, 97-140 already implement
  6.603 / 6.605), groundflow_phreatic_level_n (GROUNDFLOW_PHREATICLEVEL_N, groundfl.cc:522); groundflow_phreatic_bounda
  is not a pure alias (GNU -method1/-method2 when present, Pro switch default -yes). force_volume_element_group_0 is a
  PDF header artifact of 6.504 (force_volume_element_group already translated, database.cc:10464). slide.cc has no
  damping / user hook; groundfl.cc:653-690 for _project / _only; force.cc:226/265-278 for force_edge_diagram /
  force_gravity_geometry; general.cc:165-215 + top.cc:1734-1746 element_dof_initial* to extend by group;
  group.cc:140-153 / elem.cc:208-211 for group_time_fill; truss.cc:74 young; beam.cc has no plasticity;
  conspr.cc:170-185 direction automatic; strain_volume_absolute_time / _element registered (database.cc:8666) but
  never applied (only post.cc:164 reads them; Pro example other/strain_volume1); bounda_time_smc is a PENDING
  error (bounda.cc:458-466).
- 2026-11-13 U1 (route inline): three db_number aliases. Professional probe: the same input stores the same
  values under the manual names (it needs group_type -contact_spring, post_calcul -groundflow_pressure
  -total_pressure and control_print_database -separate_index; GNU accepts the first two).

- L-T5-bounda (2026-11-14) Owner, on "implement the switch, decide the default after measuring the
  Professional": "sí, me parece bien". Measured (column of three quad4, level y = 1.5): Professional absent
  and -yes give the same .dbs (nodes on/above the level bounded, head rho*g*z_L); -no bounds nothing. GNU
  without the record already gave the same total pressures (groundfl.cc free-surface block), so the
  default stays -yes and needs no owner decision. With bottom bounda_dof (Professional head -2.5 = GNU
  dynamic -1 + rho*g*z_L): -yes equal heads and totals; -no equal heads, total pressure at the node just
  above the level -0.5 GNU vs 0 Professional (post only): DIVERGENCES PHREATIC-BOUNDA-NO-POST.
  Note: awk under the es_ES locale reads "-0.5" as -0; compare .dbs numbers with LC_ALL=C.

- L-T5-fveg0 (2026-11-14) Owner, on discarding force_volume_element_group_0 as a PDF artifact:
  "decarta el nombre de la lista".


- L-T5-only (2026-11-15) Owner: "Mide lo qeu hace pro y hablamos"; after the measurement (Professional removes
  -groundflow from the group_type of the groups not listed in groundflow_phreatic_level_multiple_element_group and
  deletes their group_groundflow_* records; a group left without a type stops with "Error detected for data item :
  group_type") and the proposal (same database edit at start, clear error message), the owner said "si".
  Coupled column: database edit and total pressures equal; displacements already differ without the record
  (pre-existing coupled-mechanics difference, not this unit).
- L-T5-project (2026-11-15) Owner: "si" to measuring groundflow_phreatic_project on the Professional first, then
  "si, claro" to fixing the inclined single level first and adding -project on top. Measured (2x3 quad4, level
  0 1 2 2.5): Professional bounds the dry nodes to head rho*g*z_L(x) and water seeps; GNU bounded pres 0 (no flow,
  pre-existing, not recorded). -project -yes: head of the nearest point of the level line, clamped to its ends.
  Fix: pres dof = head - rho*g*z_ref (z_ref first table height), same totals as the Professional.
- L-T5-multiple (2026-11-16) Owner: "si, mide" (measure the inclined multiple level on the Professional), then
  "si, aplícalo". Measured: Professional gives the single-level seepage result; GNU bounded pres 0 (no flow).

- L-T5-slide (2026-11-17) Owner: "¿y no podemos hacer amortiguamiento viscoso y documentar la divergencia? es que
  sencillamente lo de Pro no es lo que dice que es. Puedo simular ese "amortiguamiento" con un resorte y si lo
  implementamos tal cual nos quedamos sin amortiguamiento real". Measured: Professional force (k + c/dt)*u (u total
  displacement), dt 0.01/0.02; control -no = undamped; slide_user -yes stops "Error: routine
  user_slide_friction_force is not programmed.". Implemented viscous -c*v.

- L-T5-dofinit (2026-11-17) Owner: "si, y cuando vuelvas dime exactamente en qué consiste la decisión", then
  "B documentado". Measured (bar, node_dof 2 = (1+past)/2): Professional order element full > group full > element
  specific > group specific. Option B: element full > element specific > group full > group specific (differs
  only when an element specific record meets a group full record).
- 2026-11-18 L-T5-diagram (owner "si" to measuring U5 in the Professional first). force_edge_diagram measured
  (/tmp/fd.OKmy, quad4 top edge load): factor per edge node from that node's vertical displacement of the previous
  iteration (one step = 2 iterations: -0.09; two columns E 10/100 converge to -0.1/1.1 and -0.01/1.01, so per node,
  not per integration point); 0 outside the diagram. No divergence: GNU area.cc FORCE_ELEMENT_EDGE node_factor *=
  table_xy(diagram, new_dof vertical dis). RED "I don't know what to do with : force_edge_diagram"; GREEN
  force_edge_diagram1/2 (Professional values to 12 digits), _error1 (no materi_displacement). Route inline.
- 2026-11-18 L-T5-gravgeo: force_gravity_geometry measured (/tmp/fg.dsrG): an element gets the gravity only when all
  its nodes are inside (quad4 sweep x = 0.9 none, 1.0001..1.9 first column, 2.5 both; trusses the same). A first sweep
  with the quadrilateral points in ring order was a crossed quadrilateral (manual order is the quad4 one). GNU
  force_gravity_element_inside() in materi.cc / truss.cc. RED: the keyword was in the 58 pending names of the gap probe
  (not re-run on a base build). Route inline.
- 2026-11-19 L-T5-planes: user "si" (measure U6 in the Professional first). Probes /tmp/csp.5xkX (conspr7 + horizontal
  load 0.5 per top node; sheared top block, centroid direction (-0.447,-0.894)). Professional: planes -no -yes -yes ->
  (0,-1), -yes -no -yes -> (-1,0), (0,-1) with switch_y -no -> "Error, please specify a valid contact spring direction";
  three switches required also in 2D. While comparing, three pre-existing GNU divergences showed up and were fixed:
  f = 0 without a friction record (manual 6.604: very large; Professional tangential -0.5, disx 8.5 vs GNU 1.5e12),
  2D first tangent along z (Professional t1 = (-n_y, n_x)), scalar_dabs in the spring matrix (inclined split
  0.1266/2.1094 vs GNU 0/2.236). All five cases now match to about 11 digits. Not changed: without any direction
  record the Professional stops, GNU takes the automatic one.
- 2026-11-19 L-T5-trussdiag: Professional (/tmp/trd.ATEl, length 2, diagram -0.01 -20 0 0 0.01 10 0.02 15): F follows
  the diagram of the total strain (5/12/15/-10, 0 outside), the same after two steps, after unloading and with
  materi_velocity_integrated; young and initial force ignored; an odd value count gives force 0 silently (GNU stops,
  input check as in force_edge_diagram).
- 2026-11-19 L-T5-conspr-owner: user "déjalas así, no?" -> the three contact spring fixes (unlimited friction
  without a record, 2D in-plane tangent, signed matrix) stay; independent verifier at the end of T5.
- 2026-11-19 L-T5-timefill: user "si" (U7, measure first). Professional (/tmp/tf.hgFF, two stacked quad4, upper
  group_time_fill 0.5 2.5 100): density factor 0/0.25/0.5/0.75/1 at t = 0.5..2.5 (end of step), displacement since
  birth; Mohr-Coulomb elastic while filling, plastic from t = 2.5; static phreatic level y = 2: the whole response of the
  filled element scales with f (pore pressure too). GNU: group_time_fill_factor() in group.cc, used by
  group_time_alive, get_materi_density, the materi.cc pressure factor and plasti_rule. The plastic tail differs from the
  Professional by the same amount as a plain group_time (pre-existing Mohr-Coulomb difference, not U7).
  Incident: during the U6 suite 456 repo files on the NAS turned mode 000 (content intact, git clean); restored with
  chmod 644, cause unknown (not the repo scripts).
- L-T5-strainvol (2026-11-20) U8 Professional probes (/tmp/sv.Z6ON, hex8 E 1e6 nu 0.2 total_linear, rollers):
  kernel = eigenstrain e0 = eps_v/3 on xx, yy, zz (zz also in plane strain), total value of the time table at the
  current time (back to 0 -> displacement 0), stress = D(eps - e0 I). Fully fixed hex, 3 %: sxx = syy = szz =
  -16666.67 = -E/(1-2nu)*0.01. Relative: eps_v = pct/100. Absolute: eps_v = dV / sum of the initial volumes of the
  selected set (no selector, 0.6 on V 12 -> uniform 5 %, post absolute 0.6). Free 3D: post absolute 0.1234 exact
  (strain_volume1). Post absolute = sum of tr(eps)*V (small strain), relative = 100*absolute/initial.
  Selection: no selector -> all elements; element_group; geometry = all nodes inside (centroid alone not enough).
  Pro deviations from the manual: strain_volume_element 10 1 2 uses ONLY the first element (2 1 -> initial 8 =
  element 2); relative and absolute both given -> absolute wins silently (manual: "Exactly one ... not both");
  the same element in two indices -> last index wins (6 % not 9 %), post of both indices shows that element.
- L-T5-strainvol-owner (2026-11-22) owner, verbatim: "si, perfecto" to "¿Implemento U8 siguiendo el manual en los casos 1
  y 2, con divergencias documentadas, y como Professional en el 3?". GNU: every listed element, both tables = error,
  last index wins. Kernel in materi.cc (strain_volume_setup / _in_set / _eigen): the strain passed to set_stress
  loses eps_v/3 on the diagonal; ELEMENT_STRAIN_VOLUME = integral of tr F - 3 (preallocated in top.cc for the
  parallel loop); post.cc sums it per set. All 10 tests match the Professional to 1e-8 (plane strain 0.0454000102547,
  geometry 0.100124059301, two indices 0.139875940699 / 0.460124059301).
- L-T5-baseline (2026-11-22) Professional 02-08-2026 probes (1D truss, bounda_dof 10 1 -accx, bounda_time 10 0 1 .5 1
  1 1 2 0): parameters c0 c1 c2 in tau = (t-ts)/(te-ts), record bounda_baseline_correction_parameters <bounda index>;
  given parameters are used as is. The Professional adds p(tau_k) to the bounda_time points with ts < t_k <= te
  (first point not corrected): with (-1 0 0) the applied acceleration is 1-2t on [0,.5], velocity drifts to 0.2.
  Its fit is dt independent but not reproduced ([0,2]: -0.8182607915309 -1.375947429918 2.186213875160; a = t gives
  -2 0 0). GNU without correction already equal (velocities 0.5 / 1 / 1.45).
- L-T5-baseline-owner (2026-11-22) owner, verbatim: "si" to "¿Implemento la corrección según el manual (p(τ) sumada a
  la señal en todo [ts,te], coeficientes por mínimos cuadrados de las velocidades de los pasos), manteniendo la
  convención de coeficientes de Professional y documentando la divergencia?". GNU: bounda.cc bounda_baseline_load,
  fit once with the step size of the first step in the interval (backward v += a*dt), only -acc* dofs corrected.
  [0,2] gives -0.9027362126125 -1.055204571298 2.050629834893 (sum v^2 21.68 -> 6.34e-4).
- L-T5-beamplasti (2026-11-22) Professional 02-08-2026 probes (3D cantilevers E 1 Izz 4): the limits only clip the
  components of element_beam_force_moment, stiffness stays elastic: tip disy 0.1 with mz1 limit 0.6 gives mz1 -0.6
  and fy 1.2 (out of equilibrium); force control F = 0.1t with mz2 limit 0.6: m1 -3.2 instead of -2 at F 1, then
  divergence (rotations -9.6, m1 38.4).
- L-T5-beamplasti-owner (2026-11-22) owner, verbatim: "pues prefiero un modelo plástico real. Lo que hay ahora es,
  por decirlo suavemente, una porquería. Bastaría con un método iterativo, en el que la barra se subdividiese en el
  punto en el que el ciclo de carga dé una rótula. Creo que eso no será difícil de implementar (obliga a un
  remallado al introducir la rótula) y nos daría soluciones más creíbles". After the parent noted the limits are per
  element node (hinges at element ends, no remesh): "tienes razón, había pensado en la formación de rótulas
  intermedias. No es necsario el remallado. Lo normal es mallar con malla densa donde se esperen rótulas. En este
  caso, si la densidad es suficiente, podríamos pensar mejor en mover ligeramente los nodos para que se pongan sobre
  la rótula, siempre que el movimiento sea pequeño y no distorsione la malla"; parent recommended no node relocation
  (nearest node error <= q h^2/8, history of the moved elements), owner: "si".
  GNU: beam.cc beam_plastic_hinges / beam_plastic_plane (closest point projection on |m1|<=Y1, |m2|<=Y2,
  |m1+m2|<=L V). Tangent condensed only for limits reached at the previous converged step: condensing the hinges
  of the first iterate gave a false mechanism (Lapack singular, predictor overshoot), condensing those of the
  previous iterate made Newton cycle (period 4). Hinge-forming step converges linearly (20 iterations: 5e-4 in
  deflection), exact after the next step.
- L-T5-verify (2026-11-22) end-of-T5 independent verifier, baseline 9777b2d, candidate 74ce053: U1..U10 PASS (committed
  tests plus own probes). Differential of 997 test-2014/test-2026 inputs: exit codes identical; one blocker, groundfl.cc
  single and multiple phreatic bounds wrote -0 (negative gravity times a zero height difference) where the baseline
  wrote 0 (post_calcul_static_pressure1, post_calcul_safety_default1); fixed with + 0., both .dbs identical to the
  baseline again, groundflow_phreatic* 10/10, suite 1301 OK / 0 FALLO. Advisories: force_element_edge_water stores
  the -yes code in a double (shifts with every new name, pre-existing); post_strain_volume_* records in
  post_count_data1 (U8 output); contact springs have almost no suite coverage outside the U6 tests.
- L-T2-close (2026-11-22): owner "si" to closing T2. C110 had already been resolved (L-C110, 823d58b), so the only
  open item was the 6.951 default. Pro probe (elastic velocity + one iteration with temperature-dependent conductivity):
  ratio 2.56462e-15 without the record = -materi_velocity, 0.980198 with -condif_temperature; GNU without the record
  0.980198 (all principal dofs). Fix inline in post_node_rhside_fixed_free: no record + materi_velocity ->
  {-MATERI_VELOCITY}. RED 0.980198, GREEN 6.5e-13. Tests post_node_rhside_ratio_dof_type1/2. T2 closed.
- L-MC (2026-11-23) owner: "si, adelante" (MC plastic difference follow-up), then "mantenemos GNU y documentamos".
  Measurement: one quad4, homogeneous plane strain, total_linear, E 1000, nu 0.25, phi 30 deg, c 1, final strain
  exx -0.04 eyy 0.01 exy 0.01 at constant velocity, 1/4/16/64 steps. GNU = closed-form principal face return to all
  digits, step independent (phiflow 0: sigxx -35.9458 sigxy 4.77832 sigzz -12 eppzz 0; phiflow = phi: sigxx -44.4136
  sigzz -14.8915). Pro: 1 step same principal stresses, rotated axes (sigxy 4.83215); refining drifts away (64 steps
  sigzz -13.3526, eppzz 0.00159 on a face). Documented as MC-RETURN-PRO (DIVERGENCES, SEGUIMIENTO Descarte, user and
  developer pages; K15 note). Guard tests mc_return_exact1/2 (GREEN by design, the Pro fails both). No source change.

