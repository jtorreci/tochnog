# control_safety_slip

## Description

`control_safety_slip index switch` (manual 6.369): slip-surface safety
factor calculation in the time steps of control index `index`. The slip
surfaces are given with `safety_slip_multi_linear` (6.1034-6.1037) and
`safety_slip_combined_linear` (6.1021-6.1024, same records with
`combined_linear`; the data are independent pieces x0 y0 x1 y1, x2 y2 x3
y3, ... instead of a polyline), and the circles
`safety_slip_circle_grid_*` (6.1009-6.1014) and `safety_slip_circle_line_*`
(6.1015-6.1020):

```
safety_slip_circle_grid_middle index x_first y_first [x_last y_last]
safety_slip_circle_grid_middle_n index n     (n x n middles on the rectangle)
safety_slip_circle_grid_radius index r_first [r_last]
safety_slip_circle_grid_radius_n index n
safety_slip_circle_grid_segment_n index n    (points per circle, default 90)
safety_slip_circle_grid_result index x y r safety_factor   (output)
```

`safety_slip_circle_line_*` is the same with n middles on the line from
the first to the last middle.

```
safety_slip_multi_linear index x_first,0 y_first,0 x_first,1 y_first,1 ... [x_last,0 y_last,0 ...]
safety_slip_multi_linear_n index n          (number of lines between first and last)
safety_slip_multi_linear_segment_n index n  (points per line piece, default 10)
safety_slip_multi_linear_result index x_0 y_0 x_1 y_1 ... safety_factor   (output)
```

In 3D (`number_of_space_dimensions 3`) the slip surfaces are ellipsoids
(`safety_slip_ellipsoide*`, 6.1025-6.1029) and Surfer grids
(`safety_slip_grd*`, 6.1030-6.1033):

```
safety_slip_ellipsoide index xm ym zm b1x b1y b1z b2x b2y b2z a b c [12 values of the last ellipsoid]
safety_slip_ellipsoide_n index m          (m values of every parameter with first != last)
safety_slip_ellipsoide_segment_n index n  (n x n angles, default 90)
safety_slip_ellipsoide_method index -safety_slip_ellipsoide | -materi_displacement | -materi_velocity | -materi_velocity_integrated
safety_slip_ellipsoide_result index xm ym zm e1x e1y e1z e2x e2y e2z a b c safety_factor   (output)

safety_slip_grd index -yes                (reads the Surfer ascii grid <index>.grd)
safety_slip_grd_segment_n index n         (n x n points per grid part, default 10)
safety_slip_grd_method index -safety_slip_grd_direction | -materi_displacement | -materi_velocity | -materi_velocity_integrated
safety_slip_grd_method_direction index dx dy dz
safety_slip_grd_result index safety_factor   (output)
```

## Effect on the calculation

- At the end of every time step of the control index (switch `-yes`),
  for every slip surface, sampled in points k with weight A_k:
  Fs = sum A_k tau_mc_k / sum A_k |tau_k|, with
  tau_mc = c - sigma_n tan(phi) from `group_materi_plasti_mohr_coul`
  (or `_mohr_coul_direct`) of the element group, sigma_n = n.sigma.n and
  tau the shear stress along the surface.
- The stress in a point is the average over the elements containing the
  point of the element nodal stresses (`element_dof`) interpolated in the
  element.
- Points outside the mesh, in elements without a Mohr-Coulomb record, or
  with sigma_n > 0 (tension) do not count. A surface without counting
  points gets no result.
- Circles: `segment_n` points at the angles 0, 360/n, ...
  (counterclockwise from +x), outward normal, weight 2 pi r / n; the
  surfaces of an index are ordered radius first (outer loop), then the
  middles (x fastest in the grid). The parts of a circle outside the mesh
  do not count.
- Line pieces: `segment_n` equidistant points per piece (the pieces of a
  polyline do not share points), trapezoid weights.
- With `safety_slip_multi_linear_n` the record holds the first and the
  last polyline (same number of points); n polylines are interpolated
  equidistantly between them. Without it the record is one polyline.
- `safety_slip_multi_linear_result` gets the critical polyline of the
  index and its safety factor (`safety_slip_combined_linear_result` the
  critical pieces); every piece is weighted with its own length.
- The evaluated surfaces are written as `safety_slip_surface_areas`,
  `_elements`, `_nodes`, `_critical` (the critical surface of every
  index of every type),
  `_factor`, `_stresses` (sigma_n |tau| per point), `_stress_limit`
  (tau_mc per point), `_data_index`, `_data_name` (index = surface number;
  the types are numbered in the Professional order: in 2D circle_grid,
  circle_line, combined_linear, multi_linear; in 3D grd, ellipsoide). The
  line and circle types are 2D surfaces and are not evaluated in 3D; the
  ellipsoide and grd types need 3D.
- Ellipsoid: e1 = base1 normalized, e2 = base2 orthogonalized to e1
  (Gram-Schmidt) and normalized, e3 = e1 x e2 (the result record writes
  e1 and e2). n values theta_i = 2 pi i/n and phi_j = 2 pi j/n; point
  i*n+j at r = (a sin theta cos phi, b sin theta sin phi, c cos theta) in
  that base, weight |r_theta x r_phi| (2 pi/n)^2, normal the normalized
  gradient (r1/a^2, r2/b^2, r3/c^2). The surface is covered twice (theta
  over 2 pi), which cancels in Fs. Points with zero area (the poles) get
  sigma_n 0, tau 0 and limit c. Slip direction: e1 projected on the
  tangent plane (`-safety_slip_ellipsoide`, default), or the last
  calculated displacement (`-materi_displacement`), the velocity
  (`-materi_velocity`) or `materi_velocity_integrated` at the point
  projected on the tangent plane; tau = traction . slip direction. With
  `safety_slip_ellipsoide_n m` every one of the 12 parameters with first
  != last varies over m values, all combinations (the first varied
  parameter is the outer loop): m^k surfaces.
- Grd: the file `<index>.grd` in the working directory (DSAA; nx ny; xmin
  xmax; ymin ymax; zmin zmax; then z for all x of the first y, of the
  second y, ...) defines z(x,y) on the nx x ny grid. Every part between
  four grid points gets n x n points (bilinear, the parts do not share
  points); its normal is the normalized sum of the cross products of its
  triangles 00-10-01 and 10-11-01, its area (the sum of the two triangle
  areas) is distributed with 2D trapezoid weights. Slip direction:
  `safety_slip_grd_method_direction` projected on the tangent plane
  (`-safety_slip_grd_direction`, default; required then) or a dof as for
  the ellipsoid. A missing or malformed file stops with an error.
- 3D surfaces are stored as triangles in `safety_slip_surface_elements`
  (3 point numbers each); 2D surfaces as point pairs followed by -none.
- Not evaluated in control indices without the record; switch `-no` turns
  it off.
- Switch = a number N (1, 2, ...): N passes. Every pass evaluates all
  surfaces; after a pass the zone of every index (first..last line, or
  middles and radii) is replaced by the critical surface so far plus and
  minus half a step (after the first pass) or a quarter of a step (after
  the later passes) of every parameter, with the same number of
  variations (as the Professional 02-08-2026); for an ellipsoid every
  varied parameter. A grd surface is evaluated again in every pass. The
  surfaces of all passes are kept.
- `safety_slip_set index index_0 index_1 ...` (6.1038): the indices of the
  slip surface records (of all types) that form a set;
  `safety_slip_set_result index safety_factor` (6.1039) gets the minimum
  over the set. With sets, only the critical surface of every set is
  flagged in `safety_slip_surface_critical`.
- GiD (`control_print_gid`, 6.310): the slip surfaces are written as an
  extra mesh `SafetySlip` (line elements of 2 nodes between the surface
  points; in 3D triangles of 3 nodes on the ellipsoid and grd surfaces,
  mesh `SafetySlipTriangles` when line surfaces are also present; node and
  element numbers after the highest node and element)
  with the nodal results `safety_factor` (of the surface),
  `safety_factor_local` (tau_mc / |tau|), `safety_sign` (sigma_n) and
  `safety_sigt` (|tau|). `control_print_gid_safety_slip_critical index
  -yes` plots only the critical surfaces. As every GiD mesh, the slip mesh
  is written at the first GiD print (with `control_print_gid -separate`
  at every print); later steps write the results of the same surfaces.

Same results as the Professional 02-08-2026 in elastic states (to the
solver accuracy), except `safety_slip_combined_linear_n`: the
Professional weights all pieces with the first piece (DIVERGENCES
SAFETY-SLIP-COMBINED-WEIGHTS) and sets holding index 0 (the Professional
writes no result for them and can flag a surface that is not the minimum
of the set: SAFETY-SLIP-SET-INDEX-0), `safety_slip_ellipsoide_n` (the
Professional weights every surface with the points of the first one:
SAFETY-SLIP-ELLIPSOIDE-WEIGHTS) and `safety_slip_grd_method
-materi_displacement` (refused by the Professional:
SAFETY-SLIP-GRD-METHOD). The 3D dof slip directions are interpolated by
the GNU point search to about 1e-7.

## Example

```
control_safety_slip 10 -yes
safety_slip_multi_linear 0 0.5 2. 3.5 2.
```

Uniform sigxx = sigyy = -1, sigxy = 0.5, c 1, phi 30 deg:
`safety_slip_multi_linear_result 0 0.5 2 3.5 2 3.1547005384`.

Self-checking inputs in `validation-suite/test-2026`:
`safety_slip_multi_linear1` (analytic Fs, weights),
`safety_slip_multi_linear2` (points outside the mesh, polyline),
`safety_slip_multi_linear3` (tension points),
`safety_slip_multi_linear4` (_n, _segment_n, critical line),
`safety_slip_time1` (evaluated every step),
`safety_slip_error1` (first/last lines of different size),
`safety_slip_combined_linear1` (pieces, weights),
`safety_slip_combined_linear2` (_n, numbering over the types),
`safety_slip_circle_grid1` (full and half circle),
`safety_slip_circle_grid2` (grid of middles and radii, critical),
`safety_slip_circle_line1` (line of middles, numbering),
`safety_slip_circle_error1` (middle without radius),
`safety_slip_set1` (sets over the types),
`safety_slip_refine1` (3 passes of a line zone),
`safety_slip_refine2` (2 passes of a circle line),
`safety_slip_error2` (switch 0),
`safety_slip_gid1` / `safety_slip_gid2` (GiD mesh and results, all or
critical surfaces; checked by build_safe),
`safety_slip_ellipsoide1` (sphere, weights, poles, triangles),
`safety_slip_ellipsoide2` (rotated ellipsoid, normalized bases),
`safety_slip_ellipsoide3` (Gram-Schmidt base),
`safety_slip_ellipsoide4` (_n tensor product, critical surface),
`safety_slip_ellipsoide5` (-materi_displacement, -materi_velocity),
`safety_slip_ellipsoide6` (2 passes of an ellipsoid zone),
`safety_slip_grd1` (flat grid n 3, GiD triangle mesh checked by
build_safe), `safety_slip_grd2` (inclined grid n 2, point layout),
`safety_slip_grd3` (twisted grid n 3 and n 2),
`safety_slip_grd4` (grd -materi_displacement),
`safety_slip_3d1` (numbering of the 3D types),
`safety_slip_3d_error1` .. `4` (missing method direction, missing or
malformed grd file, parallel bases). The Surfer grids `0.grd` .. `4.grd`
lie next to the inputs.
