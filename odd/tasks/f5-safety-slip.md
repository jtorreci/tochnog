# F5 T6 safety_slip: slip-surface safety factor (control_safety_slip)
Objective: implement the Professional slip-surface safety factor (manual 6.369, 6.310, 6.1009-6.1039) in the GNU,
2D first, 3D (ellipsoide, grd) as a second stage. Branch dev/f5-keywords; delivery: work-unit commits (ask-on-risk);
runner scripts/build_safe.sh (validation-suite/test-2026 self-checking inputs, Professional 02-08-2026 as oracle).

## Specs
- S1 Scope (owner 2026-11-22, verbatim "si" to "¿Te parece bien ese alcance, con el 2D primero y el 3D como segunda
  etapa?"): stage 1 = 2D: the four surface types (safety_slip_circle_grid_*, safety_slip_circle_line_*,
  safety_slip_multi_linear*, safety_slip_combined_linear*), safety_slip_set / _set_result, the automatic refinement
  (control_safety_slip N), the *_result records, the internal safety_slip_surface_* records and the GiD output
  (control_print_gid_safety_slip_critical). Stage 2 = 3D: safety_slip_ellipsoide*, safety_slip_grd*, measured first.
- S2 Safety factor (manual 6.369): "Fs = ∫ τmc dA / ∫ τ dA, where τmc is the maximum possible shear stress according
  to the mohr-coulomb condition using the actual normal stress, τ is the actual shear stress". Measured: per sample
  point tau_mc = c - sigma_n tan(phi), |tau|; points with sigma_n > 0 and points outside the mesh do not count; no
  counting point = no result record. Available for group_materi_plasti_mohr_coul, _mohr_coul_direct, _druck_prag,
  _hypo_wolffersdorff.
- S3 Stresses at a sample point: average over the elements containing the point of the element nodal stresses
  interpolated in that element (Professional, exact to 1e-12).
- S4 Discretization: line pieces segment_n points (default 10) with trapezoid weights, pieces do not share points;
  multi_linear = polyline, combined_linear = independent point pairs; circles segment_n points (default 90) at
  theta = 0, 360/n, ... counterclockwise from +x, weight 2 pi r / n. first/last data and _n: equidistant
  interpolation; circle grid n x n middles (x fastest) times radii (outer loop); circle line n middles.
- S5 control_safety_slip index switch: -yes = one evaluation per time step of that control index; a number N = N
  passes, each new pass with middles and radii in critical +/- half step, all surfaces kept. Not evaluated in a
  control index without the record.
- S6 Results: <type>_result index = critical surface of that index + Fs; safety_slip_set index i0 i1 ...,
  safety_slip_set_result index Fs = minimum over the set; critical flag: one global critical without sets, one per
  set with sets. Corrected by L-T6-U2: without sets the critical surface of every index of every type is flagged;
  surfaces are numbered over the types in alphabetical order (circle_grid, circle_line, combined_linear, multi_linear).
- S7 GiD (6.310): "If switch is set to -yes, then for a safety analysis with control_safety_slip only the critical
  slip surface will be plotted. Default, if switch is not set, all evaluated slip surfaces will be plotted ...
  Furthermore, always the normal stresses and shear stresses on the slip surfaces will be plotted." Professional
  writes mesh "safety_slip_surface_nodes" (Linear, 2 nodes) with results safety_factor, safety_factor_local,
  safety_sign, safety_sigt.

- S8 Ellipsoide (manual 6.1025-6.1029, measured L-T6-U6-measure): record 12 values (middle x y z, base1 x y z,
  base2 x y z, a b c) for the first ellipsoid, optionally 12 more for the last; e1 = base1 normalized, e2 = base2
  orthogonalized to e1 (Gram-Schmidt) and normalized, e3 = e1 x e2; the result record writes the normalized bases.
  Points: segment_n n (default 90) values of theta_i = 2 pi i/n and of phi_j = 2 pi j/n (point k = i*n + j), local
  r = (a sin t cos p, b sin t sin p, c cos t), global = middle + r1 e1 + r2 e2 + r3 e3 (the surface is covered twice,
  which cancels in Fs); weight |r_t x r_p| (2 pi/n)^2; normal = normalized gradient (r1/a^2, r2/b^2, r3/c^2) in the
  e basis; slip direction (method -safety_slip_ellipsoide, default) = e1 projected on the tangent plane, normalized;
  tau = traction . slip direction; zero-area points (poles) get sigma_n 0, tau 0, limit c. Methods
  -materi_displacement / -materi_velocity_integrated (last calculated displacement at the point) and -materi_velocity
  (velocity), projected on the tangent plane. safety_slip_ellipsoide_n m: each of the 12 parameters with first !=
  last varies over m values, all combinations (first varied parameter outer loop, last varied parameter inner loop);
  m^k surfaces for k varied parameters. Result record: the 12 parameters (normalized bases) + Fs. Elements record:
  triangles (i,j),(i,j+1),(i+1,j) and (i,j+1),(i+1,j+1),(i+1,j), indices wrapping in i and j, 3 values each without
  -none. Professional values (3x3x3-node hex8 cube [0,4]^3, uniform sxx = syy = szz = -1, sxy = 0.5, c 1, phi 30 deg):
  sphere (2,2,2) r 1 n 8: 10.81644604373; ellipsoid (2,2,2) base1 (1,1,0) base2 (-1,1,0.5) a b c 1.5 1 0.8 n 8:
  9.298529656311; base1 (1,0.2,0) base2 (0,1,0) sphere: 8.366977774663; sphere with -materi_displacement: 9.032320731610;
  first (2,2,2) r 1 / last middle x 2.5 and a 1.2, _n 3: 9 surfaces 10.816446, 9.771643, 9.013850 (repeated per
  middle x), critical middle x 2.5 a 1.2.
- S9 Grd (manual 6.1030-6.1033, 3D only): safety_slip_grd index -yes reads the Surfer file <index>.grd (DSAA, nx ny,
  xmin xmax, ymin ymax, zmin zmax, then z for all x of the first y, then the second y, ...) from the working
  directory: a surface z(x,y) on the nx x ny grid, (nx-1)(ny-1) parts. Every part gets segment_n x segment_n points
  (default 10, bilinear interpolation of the part corners, parts do not share points); the part normal is the
  normalized sum of the normals of its two triangles (corner 00-10-01 and 10-11-01), constant for all points of the
  part; the part area (sum of the two triangle areas) is distributed with 2D trapezoid weights (n 3: corner 1/16,
  edge 1/8, middle 1/4 of the area). Slip direction: safety_slip_grd_method -safety_slip_grd_direction (default)
  takes safety_slip_grd_method_direction projected on the tangent plane; the Professional refuses the method
  -materi_displacement (db error on safety_slip_grd_method) and stops with a db error on
  safety_slip_grd_method_direction when it is missing. Result: safety_slip_grd_result index Fs. Professional values
  (same cube, uniform sxx = syy = szz = -1, sxz = 0.5): 0.grd flat z 2 over [0.5,3.5]^2 n 3: 3.154700538384; 1.grd
  3x2 (x 0.5..3.5, y 1..3, z rows 1.5 2 2.5 / 2 2.2 2.4) n 2: 3.821906810681 (point weights 0.7864226581436 and
  0.7758425062193, sigma_n -1.21504 in part 1); 2.grd 2x2 (x 0.5..2, y 1..3, z 1.5 2 / 2 2.2) n 3: 3.840082933030.
- S10 Numbering and GiD in 3D: types in alphabetical order circle_grid, circle_line, combined_linear, ellipsoide,
  grd, multi_linear (assumption, to be confirmed by a mixed probe); GiD 3D slip mesh of triangles (ElemType
  Triangle, Nnode 3). Corrected by L-T6-U6-order: in 3D the Professional numbers grd before ellipsoide and skips
  the line types; GNU: 2D types only in 2D, grd then ellipsoide in 3D.

## Tasks
- [x] U1 control_safety_slip -yes + safety_slip_multi_linear (+ _n, _segment_n, _result) + point stresses + Fs +
  safety_slip_surface_* records. Route: inline (dependent units, no writer reason; map L-T6-map). RED: base "I do not
  know : control_safety_slip". GREEN safety_slip_multi_linear1..4, safety_slip_time1, safety_slip_error1 (analytic,
  = Professional in the elastic probes to solver accuracy). Commit 0964a37; suite 1307 OK / 0 FALLO.
- [x] U2 safety_slip_combined_linear (+ _n, _segment_n, _result), critical per index, numbering over the types.
  Route: inline. RED: base "I do not know : safety_slip_combined_linear". GREEN safety_slip_combined_linear1/2 (analytic);
  divergence SAFETY-SLIP-COMBINED-WEIGHTS.
- [x] U3 safety_slip_circle_grid_* and safety_slip_circle_line_* (middles, radii, _n, _segment_n, _result).
  Route: inline. RED: base "I do not know : safety_slip_circle_grid_middle". GREEN safety_slip_circle_grid1/2,
  safety_slip_circle_line1 (analytic; the Professional gives the same values on these inputs: 4.954385746652,
  4.893695691055, 3.872934051177), safety_slip_circle_error1. U2 and U3 in one commit (same files); suite 1313 OK / 0 FALLO.
- [x] U4 safety_slip_set / _set_result, critical flags, control_safety_slip N refinement. Route: inline. RED: base
  "I do not know : safety_slip_set"; refinement N: db_error on base (only -yes). GREEN safety_slip_set1,
  safety_slip_refine1/2, safety_slip_error2 (analytic; Professional same values except set index 0); divergence
  SAFETY-SLIP-SET-INDEX-0. Rewrite of the driver as zones (passes over all types).
- [x] U5 GiD output of the slip surfaces, control_print_gid_safety_slip_critical. Route: inline. RED: base "I do not
  know : control_print_gid_safety_slip_critical", no SafetySlip mesh. GREEN safety_slip_gid1/2 (build_safe file
  checks: 70/63 and 30/27 nodes/elements, 4 results); divergence SAFETY-SLIP-GID-FORMAT. U4 and U5 in one commit; suite 1319 OK / 0 FALLO.
- [x] U6 (stage 2) 3D: safety_slip_ellipsoide*, safety_slip_grd* (measured: L-T6-U6-measure, spec S8-S10).
  Route: delegated writer (long-session context backstop of the parent). RED: base "I do not know :
  safety_slip_ellipsoide" / "I do not know : safety_slip_grd". GREEN safety_slip_ellipsoide1..6, safety_slip_grd1..4,
  safety_slip_3d1, safety_slip_3d_error1..4 (Professional values / analytic), GiD triangle mesh check on
  safety_slip_grd1; 2D differential (18 inputs, HEAD vs U6) identical .dbs/.msh/.res. Divergences
  SAFETY-SLIP-ELLIPSOIDE-WEIGHTS, SAFETY-SLIP-GRD-METHOD. Commit: pending (parent).

## Log
- L1 2026-11-22 owner, verbatim: "si" to "¿Arranco T6 midiendo `safety_slip` en Professional?"; then "si" to the
  scope question of S1.
- L-T6-measure (2026-11-22) Professional 02-08-2026 probes (4x4 quad4 under prescribed uniform strain, MC c phi):
  sxx = syy = -1, sxy = 0.5, c = 1, phi = 30 deg: Fs 3.154700538384 on y = 2 (= 1.57735 / 0.5); full circle r 1 in
  the uniform field 4.954385746652 (90 points); tension points excluded (piece 1 tension + piece 2 compression gives
  the compression-only value 201.1547005384 with c 100); all tension: no record; point outside: stress 0, limit 0.
  Grid 3x3 x 2 radii = 18 surfaces, surface = ir*9 + iy*3 + ix; N = 2: second zone x 2.25..2.75, y 1.75..2.25,
  r 0.75..1.25 around critical (2.5, 2, 1). Evaluated every step of a control index with control_safety_slip
  (two-step index: Fs 4.309401076768 of the last step); a following index without the record keeps the old result.
  Internal records safety_slip_surface_{areas, elements (point pairs with -none), nodes (x y z), critical, factor,
  stresses (sigma_n tau pairs), stress_limit, data_index, data_name}.
- L-T6-map (2026-11-22) read-only explorer: element nodal stresses element_dof [inol*nuknwn + stres_indx + c*nder]
  (miscel.cc element_dof_nodal_set, top.cc:1020, options_element_dof -yes default); point_el (point_el.cc:86), the
  post_point loop stops at the first element so safety_slip has its own serial loop; get_group_data (group.cc:538),
  mohr_coul phi in radians (plasti.cc), mohr_coul_direct units ambiguous, hypo_wolfersdorff phi_c in degrees;
  hook after post() in step_close (top.cc:1872); GiD print_gid_6 (print_g6.cc) writes element meshes and OnNodes
  results only; control record -yes or number pattern control_print_plotmtv.
- L-T6-U1 (2026-11-22) GNU on the Professional probes: elastic cases equal to 1e-9 (GNU solver noise in the uniform
  field, Professional exact); plastic MC cases differ because the MC states differ (p_incl, v_ml), not the safety
  evaluation. materi_stress is always required by the materi groups, so no own check.
- L-T6-U2 (2026-11-22) Professional probes: mixed types (o_mix): surfaces numbered circle_grid, circle_line,
  combined_linear, multi_linear, and every index of every type gets its own critical flag (5 flags for 5 indices),
  so the "one global critical" of L-T6-measure was wrong (it had one index). combined_linear with _n: all pieces get
  the point weights of the first piece (e_comb_n surface 0: 9.041451884343 = equal weights; correct 8.548948025,
  the Professional value without _n); GNU follows the manual (SAFETY-SLIP-COMBINED-WEIGHTS).
- L-T6-U4 (2026-11-22) Professional probes: passes are global (each pass evaluates all types and indices again,
  also the ones without variations; numbering pass-major); zone after pass 1 = critical so far +/- step/2, after
  the later passes +/- step/4 (k_ml3 3 passes 3.104571210375, k_g5 n 5, k_l4 4 passes); sets: a set holding index 0
  gives no result (k_a, k_b, k_set0), sets without 0 work (k_c, k_e); with sets the per-index critical flags are
  not written; k_set2 flagged surface 1 (circle 1, 5.21) for set_result 1 = 3.1547 (multi 1, surface 3).
- L-T6-U5 (2026-11-22) Professional GiD (s_set): mesh safety_slip_surface_nodes Linear 2 nodes, node ids after the
  mesh nodes, element ids after the elements; per line element (both Gauss points) the value of its first point:
  safety_sign = sigma_n, safety_sigt = |tau|, safety_factor_local = tau_mc/|tau|, safety_factor = surface Fs. The GNU
  writes the same values on the nodes of a SafetySlip mesh (GNU GiD 6 writes OnNodes results only).
- L-T6-U6-measure (2026-11-22) Professional probes in scratch s3/ (gen3.py: hex8 cube under uniform strain):
  el1..el5 (ellipsoid sphere, rotated ellipsoid, _n tensor product, base Gram-Schmidt, -materi_displacement),
  gr1..gr5 (flat, inclined 3x2, twisted 2x2 n 3, missing direction = db error, -materi_displacement refused).
  Values and rules in S8-S10.
- L-T6-U6 (2026-11-22) writer: all per-point records of el1..el5, gr1, gr2, gr4 equal the Professional (1e-6),
  except: (a) _n: the Professional factors use the point weights of the first surface (its own per-point stresses
  give the GNU 9.5517495368 / 8.6775003379), GNU uses each surface's weights (SAFETY-SLIP-ELLIPSOIDE-WEIGHTS);
  (b) dof directions (displacement/velocity) differ by about 1e-7 because point_el inverts hex8 iteratively, tests
  use 1e-5; dis/vel are interpolated from node_dof (dof_element_dof -no for them). The grd part normal is the
  normalized sum of the two triangle cross products (area weighted), not of the unit normals (S9 wording);
  the theta = pi pole (area 1e-16) is also treated as a zero-area point. Numbering combined_linear, ellipsoide,
  grd, multi_linear checked in safety_slip_3d1 (GNU only; the Professional order of ellipsoide/grd is an
  assumption, S10). Records safety_slip_surface_* get the length of the largest surface of the step
  (db_data_length_put; SAFETY_SLIP_MPOINT 10000 is the maximum). The .grd files are gitignored in
  validation-suite/test-2026 (only *.dat are versioned).
- L-T6-U6-order (2026-11-22) parent, Professional mixed 3D probes (s3/mx, mx2): data names circle_grid, circle_line,
  grd, ellipsoide; the multi_linear and combined_linear records give no surface in 3D; the circles in 3D lie in the
  plane z = 0 (Fs ~1e15, meaningless). GNU: line and circle types evaluated only in 2D, grd before ellipsoide in 3D
  (DIVERGENCES SAFETY-SLIP-CIRCLE-3D for the circles); safety_slip_3d1 rewritten accordingly. Parent verification
  of the writer result: diff inside the allowed surfaces, ellipsoide2 / grd2 / combined_linear2 rc 0, full suite
  1335 OK / 0 FALLO before this order change and 1335 OK / 0 FALLO after it.
