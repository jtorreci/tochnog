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
- [ ] U6 (stage 2) 3D: safety_slip_ellipsoide*, safety_slip_grd* (measure in the Professional first).

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
