# GiD family (F5 follow-up, reopened discard)

Objective: implement the GiD output records and the GiD geometry records of the Professional 02-08-2026 that the
GNU discarded at the start of F5. Branch `dev/f5-keywords`. Delivery: ask-on-risk, work-unit commits, never push.
Runner: `scripts/build_safe.sh` (`TN_SKIP_SUITE=1` for builds); Professional binary as oracle where it runs.

## Specs

- S1 Owner scope (L1, L2): "La 1 y la 2, la salida y la geometría". 1 = GiD output: the `control_print_gid_*`
  (6.294-6.315) and `print_gid_*` (6.982-6.989) names still missing. 2 = GiD geometry: the `mesh_gid_*` names
  (6.818-6.853). Not in scope: `control_print_gid_batch` (6.294). `control_mesh_gid_batch` (6.200) came into scope
  with S5 as the writer of the GiD batch file.
- S2 Context (L1): the owner works on a preprocessor that reads the GiD format.
- S3 Output records (group 1), Professional behaviour measured in L5, `print_gid_X` = global default,
  `control_print_gid_X idx` overrides it for the `control_print_gid idx` with the same index:
  `_coord`/`print_gid_coord` (coord result), `_group`/`print_gid_group` (element_group result), `_other` (boundary,
  deformation, element results), `_dof` (dof fields or -none), `_dof_calcul` (post_calcul fields or -none),
  `_element_group` (mesh blocks of the listed groups), `_node_type`/`print_gid_node_type` (-node, -node_start_refined,
  -node_deformed_mesh coordinates in the .msh), `_spring2`/`print_gid_spring2` (spring2 as 1 or 2 nodes),
  `_beam_vectors` + `_beam_vectors_normal`, `_truss_vector` + `_truss_vector_normal` (vector results),
  `_smooth_dof`/`print_gid_smooth_dof` + `_smooth_n`/`print_gid_smooth_n` (nodal smoothing, n default 10),
  `_element_mpc`, `_mesh_activate_gravity`/`print_gid_mesh_activate_gravity`, `print_gid_define`.
  Professional-suspect behaviour, decided in L7-L9: the GNU follows the manual for `_smooth_dof` with named dofs and
  for `_dof_calcul` with named fields; `print_gid_smooth_n` takes no index (as the Professional).
- S4 Geometry records (group 2): the 35 `mesh_gid_*` names of 6.818-6.853 are accepted with the manual syntax (L5).
- S5 Owner decision (L6, Q "Sin GiD instalado, ¿qué hace GNU con la geometría mesh_gid_*?"): "Escribir el .bch".
  `control_mesh_gid_batch idx -yes` writes `GidMeshBatch.bch` in the Professional format (byte-equal to the
  Professional on its 8 examples and the 2D/3D probes of L5). If `gid` is on the PATH it runs `gid -n -b
  GidMeshBatch.bch`, copies `mesh_gid.gid/mesh_gid.dat` to `<idx>.dat` and deletes the .bch and the .gid directory,
  as the Professional does; otherwise it keeps the .bch, prints a warning and the run continues.
- S6 Docs: user and developer pages, DIVERGENCES (remove DSC-GID, add any Pro-suspect rows), SEGUIMIENTO discard row
  and metric, keyword probe.

## Tasks

| ID | Specs | Route | Status | Commit |
|---|---|---|---|---|
| T0 | S1 | delegated (explorers) | done (L4, L5) | - |
| T1 | S4, S5 | delegated writer (context backstop) | done (L10) | d61231c, f8176cf |
| T2 | S3 | delegated writer (context backstop) | done (L12) | f4fe4ca |
| T3 | S3 | delegated writer | done; Q-DEF resolved (L13, 0d6b607), Q-MAG in T6 | f4fe4ca |
| T4 | S3 | delegated writer | done (L12) | f4fe4ca |
| T6 | S3 | delegated writer | done (L15) | this commit |
| T5 | S6 | inline | done (L16) | this commit |

## Log

- L1 (2026-11-24) owner: "Intentemos implementarlo. Te propondría seguir con GiD cuando hagamos estos 3. Aunque lo
  descartamos desde el principio, estoy trabajando en un proyecto de un preprocesador que admite formato GiD y nos
  puede interesar"
- L2 (2026-11-24) owner, after the three groups were listed (1 output 25 names, 2 geometry 35, 3 batch 2):
  "La 1 y la 2, la salida y la geometría"
- L3 (2026-11-24) Professional evidence so far: manual 6.200 says the mesh_gid_* records "are used to generate
  elements using GiD as preprocessor. In fact, the GiD program is called from inside Tochnog as batch". The
  distribution has partners/gid/tochnog.gid (tochnog.bas, .cnd, .mat, .bat), test/other/mesh_gid_1..8.dat and 35
  test files that use print_gid or mesh_gid records.
- L4 (2026-11-24) GNU map (explorer): GiD writer print_g6.cc (print_gid_6, ISO-8859 file), called from top.cc:1964 and
  miscel.cc:1331/1534 without icontrol (reads ICONTROL from the db, print_g6.cc:117). control_print_gid_empty, _method,
  _save_difference, _contact_spring2 are parse-only. Template: VTK options via print_option(global, control, icontrol,
  default) (miscel.cc:1269), print_vt.cc:75-111. No mesh_gid code. Docs: discard rows SEGUIMIENTO:449, DIVERGENCES
  DSC-GID (:43).
- L5 (2026-11-24) Professional probes (scratch /tmp/gidprobe.uChW): print_gid_X = global default, control_print_gid_X idx
  overrides. Measured effect of _coord, _group, _other, _dof, _dof_calcul, _element_group, _node_type, _spring2,
  _beam_vectors(+_normal), _truss_vector(+_normal), _smooth_dof/_smooth_n. No visible effect: print_gid_define,
  _element_mpc (no mpc in probe), _mesh_activate_gravity (untested). mesh_gid_* without control_mesh_gid_batch: accepted
  and ignored. With it: writes GidMeshBatch.bch (GiD batch commands for points, lines, arcs, objects, surfaces, volumes,
  sizes, structured meshing, materials GROUP_g, conditions nodes_in_*), runs `gid -n -b`, reads mesh_gid.gid/mesh_gid.dat
  (tochnog.bas output: node, element, element_group, geometry_list 1-40), copies it to <idx>.dat, deletes the .bch.
  Pro-suspect: smooth_dof works only with -all; _dof_calcul with a name suppresses it; manual 6.200 example uses
  control_mesh_gid (rejected); print_gid_smooth_n rejects the index the manual shows; volume_element_group emits
  "Data Mat Assign 3 Volume"; missing gid still exits 0.
- L6 (2026-11-24) owner, Q "Sin GiD instalado, ¿qué hace GNU con la geometría mesh_gid_*?": "Escribir el .bch
  (Recomendado)" (option text: control_mesh_gid_batch writes GidMeshBatch.bch in the Professional format; if gid is on
  the PATH it runs it and reads the mesh; otherwise it keeps the .bch and warns). Fixtures preserved in the session
  scratchpad gidfix/ (mg1..8, m2d, m3d kept.bch, fake gid script).
- L7 (2026-11-24) owner, smoothing (Pro smooths only with -all, named dofs give no change; manual: smooths the named
  dofs): "Seguir al manual (Recomendado)" -> GNU smooths the named dofs, -all smooths all; Pro-suspect divergence row
  plus guard test.
- L8 (2026-11-24) owner, _dof_calcul (Pro: a named field is suppressed like -none; manual: only the named fields are
  written, parallel to _dof): "Seguir al manual (Recomendado)" -> GNU writes only the named fields; Pro-suspect row plus
  guard test.
- L9 (2026-11-24) project decision, no owner question: `print_gid_smooth_n` takes no index (the Pro rejects the index
  the manual example shows, and every other global print_* record has none); manual erratum, not a divergence.
- L10 (2026-11-24) T1 done (writer): 36 mesh_gid_* names (the manual has 36, not 35) + control_mesh_gid_batch;
  mesh_gid.cc writes GidMeshBatch.bch byte-equal to the Professional on mg1..mg8, m2d, m3d and 6 extra probes; with a fake
  gid on PATH, 10.dat equals the Professional run and the .bch and mesh_gid.gid are removed; without gid a warning and
  rc 0. RED: "I do not know : mesh_gid_point_coord". Suite 1392 OK / 0 FALLO. Parent spot check: mesh_gid_2d1 cmp OK.
  Open: pre-existing GNU bug, control_input without unknowns fails ("Error detected for data item :
  node_start_refined", input_control_input creates node_start_refined only when nuknwn>0), so the Professional's own
  examples stop after reading 10.dat; GNU stricter than the Pro on bad assign_conditions indices; 3D circles skipped
  (the Pro writes a broken command); 2D cylinder/sphere/volume rejected at batch time. Windows makefile mirrored by
  the parent. Risk: item 6 (the real-gid branch is tested only with a fake gid).
- L11 (2026-11-24) Pre-existing GNU defect fixed (found by the T1 writer): control_input in a model without dofs did not
  create node_start_refined for the nodes read ("Error detected for data item : node_start_refined"), so the
  Professional mesh_gid examples could not read their mesh. RED observed on d61231c with control_mesh_gid_batch3, GREEN
  after the fix. Suite 1393 OK / 0 FALLO.
- L12 (2026-11-24) T2-T4 done (writer): 16 control_print_gid_* + 8 print_gid_* records, print_gid(task, icontrol),
  print_option_list() in miscel.cc, 27 tests (RED on HEAD: "I do not know", GREEN 27/27); L7 guard
  control_print_gid_smooth_dof2, L8 guard control_print_gid_dof_calcul1; smoothing, truss_vector, beam_vectors,
  -node_deformed_mesh and _element_mpc -no equal to the Professional to 1e-9 (smoothing = Jacobi pass 0.75 old + 0.25
  mean of element neighbours, bounded dofs kept). Differential HEAD vs new on 1170 inputs: 1818 flavia files identical.
  Suite 1420 OK / 0 FALLO. Open: Q-DEF the Professional writes coord and element_group by default (GNU only on request,
  to keep existing output); Q-MAG _mesh_activate_gravity parse-only (Pro default prints only the activated elements);
  _other -no removes only what the GNU writer writes (no truss_force, beam_force_moment, bounda_force, node_boundary);
  print_gid_define parse-only (no Pro effect); a spring2 solution difference GNU vs Pro noticed (not in scope).
  Risk: item 3 (print_gid signature, renumbering arrays), covered by the differential.
- L13 (2026-11-24) owner, Q-DEF ("Por defecto, el Professional escribe en la salida GiD el resultado coord y el campo
  element_group, y dibuja el spring2 con 1 nodo... ¿Cambiamos el default de GNU al del Professional?"): "Default del
  Professional (Recomendado)". print_gi.cc defaults: coord -yes, group -yes, spring2 1; every .flavia.res gains coord
  and element_group. Check control_print_gid_dof2 extended to pin the defaults.
- L14 (2026-11-24) owner, Q-MAG (Pro: by default the GiD mesh holds only the elements already active at the first GiD
  print, -yes holds all; GNU always all): "Copiar al Professional (Recomendado)". Task T6.
- L15 (2026-11-24) T6 done (writer): Professional rule measured: by default the .flavia.msh is written once at the first
  GiD print with the elements active then; each print's element results hold the elements active at that print; -yes all
  elements; -method2 active for t >= start, -method1 t > start, inactive up to time_initial. mesh.cc exports
  mesh_activate_gravity_active(). 5 tests (RED 3/5 on HEAD as expected, GREEN 5/5). Differential on 1174 inputs: only
  models with mesh_activate_gravity change (mesh_act_grav2 now an empty mesh, as the Professional writes none;
  mesh_activate_gravity_time1 now equal to the Professional). Suite 1424 OK / 0 FALLO. Pro-suspect, copied: element
  results for elements not in the mesh; empty .msh when nothing is active at the first print; method1 > vs method2 >=.
  Not in scope: the Professional writes post_point as Point elements in its GiD mesh.
- L16 (2026-11-24) T5 done: keyword probe on 8a53ace: 1137 accepted + 6 literals + 2 normalised = 1145 of 1147;
  discarded control_print_gid_batch and force_volume_element_group_0. SEGUIMIENTO metric, discard row (reopened) and 64
  GiD checklist entries; DIVERGENCES: DSC-GID narrowed to control_print_gid_batch, new GID-SMOOTH-NAMED,
  GID-DOFCALCUL-NAMED, MESH-GID-BATCH-ERRORS (class counts table not recomputed). Next: feature-end independent verifier.
