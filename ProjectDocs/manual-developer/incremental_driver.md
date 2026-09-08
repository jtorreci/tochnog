# incremental_driver — ground-truth algorithm (measured from the Professional binary 25-10-2023)

DEV-A sprint doc. Everything here is MEASURED (Pro A/B probes + objdump of
`incremental_driver*` symbols), not invented. Raw disassembly: /tmp/incdrv_dis_v2/
(annotated *_annot.txt). Pro probes: /tmp/incdrv_probes/, /tmp/incdrv_pro/.

## 1. Program-level hook (top.cc)

In the Pro, `top()` calls `calculation_start(); easy(); incremental_driver();`
(0x21fae3..ed). `incremental_driver()` (0x1a6d00):

```
if (!db_active(INCREMENTAL_DRIVER)) return;          // 1198 in Pro enum
if (!db_active(INCREMENTAL_DRIVER_COMBINED))          // auto-create: combined[i] = {i}
  for i in 0..highest_index(INCREMENTAL_DRIVER): db(PUT, COMBINED, i, {i});
file_delete("incremental_driver_result.txt");
open result file; write 2 header lines (56 + 146 chars);
precision = print_precision, scientific;
parallel_processors = processors_maximum();
parallel_sys_routine(parallel_incremental_driver);    // the whole driver
close file;
exit_tn();                                            // THE PROGRAM EXITS HERE
```

=> When `incremental_driver` is active the normal FEM loop NEVER runs. The GNU
hook must: after parsing (before the step loop), if the record is active, run the
driver, write the result file and exit with status 0.

## 2. Result file format (measured, byte-exact)

```
line1: "Results are given for each timestep on a different row:\n"   (56 chars)
line2: "time stress_xx stress_xy stress_xz stress_yy stress_yz stress_zz strain_xx strain_xy strain_xz strain_yy strain_yz strain_zz pressure void_ratio\n" (146 chars)
per combined experiment c:
  "\nresults for combined experiment %d\n"     (c = index of the combined record)
  "number of experiments %d\n"                 (length of combined[c])
  "number of timesteps %d\n"                   (ntime, per experiment)
  ntime*nexp rows: 15 scientific columns
    time sxx sxy sxz syy syz szz exx exy exz eyy eyz ezz pressure void_ratio
```
- `time` is ABSOLUTE within the combined sequence (experiments continue the clock;
  a combined experiment re-runs FROM THE INITIAL STATE — combined 1 in the syntax
  test repeats combined 0's first 1000 rows bit-for-bit).
- values format: `%.12e` separated by single spaces, trailing space, CRLF-free.
- pressure: pore pressure, compression-POSITIVE. void_ratio: only updated by hypo.

## 3. Records (Pro enum → name; GNU adds new)

1198 incremental_driver (value: -incremental_driver_oedometric/-triax/-direct_shear/
       -hydrostatic/-pure_shear), 1199 _combined, 1200 _element_group,
       1201 _experiment_force_controlled, 1202 _experiment_time,
       1204 _intergranular_strain, 1205 _isa_strain, 1207-1209 _length_x/y/z,
       1210 _void_ratio, 1211 _pressure, 1215 _stress, 1217 _volume_constant,
       1218 incremental_ntime (default 1000),
       1114-1116 group_materi_plasti_mohr_coul_direct_incremental_driver_{c,phi,psi}
       (piecewise (gamma,c/phi/psi) tables; gamma = total engineering shear strain).

## 4. Per-experiment semantics (measured)

Experiments run in the order listed by combined[c]. Timeline: experiment i spans
[experiment_time[0], experiment_time[N-1]] (absolute); prescribed value linearly
interpolated; the value is ADDITIONAL to the previous experiment's end value.

- oedometric: eps_x=eps_y=0; z: displacement-controlled compression.
- triax: x (and y if not axisymmetric): FIXED STRESS (lateral cell); z: displacement.
  For axisymmetric triax the y (theta) strain comes from the solve (measured
  undrained: eps_xx = eps_yy, i.e. equibiaxial lateral strain).
- direct_shear: y fixed, z fixed stress, x shear of top plane.
- hydrostatic: volumetric.
- force_controlled -yes: the wanted force/area replaces displacement (Newton both
  ways); all 4 targets are displacement-controlled.
- -updated_area/-updated_linear: SMALL strain (eps = delta_d/L0, lengths stay
  initial; measured: oedo eps rows are exactly linear). The oedo result is
  IDENTICAL for both memory values (290287 bytes both) — area only matters for
  force-controlled runs.
- initial equilibrium: lateral wanted stress for experiment i>0 = end value of the
  previous experiment; for the first experiment the lateral wanted force/area =
  current lateral stress (0.5*(sig_xx+sig_yy) axisym).

## 5. Material step (per iteration; measured + disassembled)

Per Newton iteration k:
1. `incremental_driver_deformation` (0x1a9d20): builds the trial deformation from
   the current iterate, computes areas (axisymmetric: 2*pi*r*L side / pi*r^2 top),
   assembles total_wanted_force_per_area per experiment type, calls
   `calculate_deformation_correction` (0x1ab240: 3x3 tangent system, LU via
   matrix_inverse, "no inverse for matrix" on singular) to correct the lengths,
   then deften -> strain/rotation via `incremental_driver_deformation_deften_decompose`
   (0x1ab8a0: C = F^T F, eigen, U = V diag(sqrt|lam|) V^T, rot = F U^-1;
   strain = symmetric part; identical math to GNU set_deften_u_rot).
2. material law:
   - mohr_coul_direct (0x1aa930): sigma_trial = sigma_old + C_matrix(young,poisson):dstrain
     (C_matrix fills the 3x3x3x3 elastic stiffness = the tangent output).
     gamma = total strain (eps_xy + eps_yx); optional tables interpolate c/phi/psi(gamma).
     beta = sin(psi)/(1-2*nu).
     Eigen-decompose sigma_trial (DESCENDING). f = mohr_coul_function_continuum:
       f = 0.5|v0-v2| + 0.5(v1+v2)*sin(phi) - c*cos(phi)     (v0>=v1>=v2)
     If f>0: gradient-flow loop (<=999 sub-iterations, each solving a 1-D
     equation for the step factor, f halved per solve, 0.1 damping, tol 1e-10,
     |sigma| tol 1e-8*array_size) that moves the principal values along a FIXED
     direction until f<=0, then matrix_from_eigenvalues.
     MEASURED direction (axisymmetric, laterals equal):
       psi=0: (1, 1, -2)/2 in principal space == mean-preserving (p-constant),
       Lode-preserving return. Verified against Pro oedo to 8 significant digits
       (sig_zz step1 -1.146480306121 vs analytic -1.146480187: Newton residual
       ~1e-7 relative; K0 = (1-sin phi)/(1+sin phi) EXACT).
       psi>0: adds a mean shift. Measured in (p,q) space (p=mean, q=sig3-sig1),
       direction (dp/dq):
         (psi=0.558, nu=0.23): -0.651016  (undrained target parameters)
         (psi=0.558, nu=0):    -0.285769
         (psi=0.3,   nu=0.23): -0.332716
       => h(beta) fit through (0,0),(0.529784,-0.285769),(0.981082,-0.651016),
       beta = sin(psi)/(1-2nu). The two target points (psi=0, psi=0.558/nu=0.23)
       are EXACT measured values.
   - hypo: void_ratio update e = max(0,(e_prev+1)*(1+evol)-1) (evol = trace-3 or
     det-1 for large deformation), then hypo_law(...) (0x182800, same signature as
     the GNU's). [DEV-C kernel dependency]
   - umat: get_group_data(group_materi_umat_parameters) + umat call.
3. groundflow (0x1ab080), per iteration, if the group is -materi -groundflow:
     capacity = GROUP_GROUNDFLOW_CAPACITY[gr] (must be > 0)
     evol  = trace(deften) - 3   (small strain)  |  det(deften) - 1 (large)
     p_new = p + evol / capacity
     water coupling: tangent diagonals += 1/capacity (2/capacity for one variant)
   MEASURED undrained invariants (target test):
     p = evol/capacity EXACTLY (capacity = group_groundflow_capacity = 1e-7)
     sigma_x_law + p = lateral applied stress (=-100) EXACTLY
     the mohr_coul_direct surface holds on the LAW stress (effective stress;
     pore pressure enters only through the lateral equilibrium + capacity law)
4. convergence: stress iterate change <= tol (measured noise 1e-8 relative in the
   rows; max_iter with "not_converged: combined experiment %d time step %d"
   messages), sigma must stay definite (new_deften diag > 0 check).

## 6. Newton structure (displacement-controlled, axisymmetric triax/oedo)

unknowns per step: lateral strain eps_l (=eps_xx=eps_yy) [triax] or fixed 0
[oedo], and p_water [undrained].
residuals: lateral: sigma_x(new) - sigma_x_wanted (+ p coupling); water: p - evol/C.
The Pro couples them in one 3x3 solve (calculate_deformation_correction) with the
material tangent (mcd's 81-double C_matrix output + water 1/capacity terms).

## 7. Verification data (Pro 25-10-2023, rc=0, deterministic — md5 stable over repeats)

- oedometric_drained_mohr_coulomb: /tmp/incdrv_pro/oedo/incremental_driver_result.txt
  (1000 rows; step1 sig_zz=-1.146480306121, sig_xx=sig_yy=-0.3526857727628;
  K0 EXACT (1-sin 0.558)/(1+sin 0.558) = 0.3076248820)
- triaxial_compression_drained: /tmp/incdrv_probes/orig_triax/ (1000 rows;
  elastic until sig_zz=-650.1426293481 then perfectly plastic at fixed
  sig_x=-200; final eps_xx=+3.784614896973e-3)
- triaxial_compression_undrained: /tmp/incdrv_probes/rep_undr/ (1000 rows;
  final p=+166.7643777647 = evol/1e-7 EXACT; sigma_xx+p=-100 EXACT;
  final law stress on the MC surface with c=10)
- incremental_driver_syntax: /tmp/incdrv_pro/syntax/ (hypo-dependent: group 0 is
  hypoplastic masin; sig_zz stays exactly -110 while shearing; combined 0 = 1000
  rows t=1e-3..1, combined 1 = 2000 rows t=1e-3..2 from the initial state)

GOTCHAS (measured):
- Pro leaves result EMPTY (0 rows, rc=0) if it fails silently mid-run — seen
  twice on first runs of triax/undr, unreproducible on repeats (transient; keep
  row-count checks when comparing).
- The Pro auto-creates incremental_driver_combined records and dumps them into
  the .dbs.
- psi/phi default from group_materi_plasti_mohr_coul_direct (phi, c, psi); the
  _incremental_driver_{c,phi,psi} tables are OPTIONAL overrides (none of the 4
  targets uses them).

## 8. ISA syntax boundary (parse-only)

`materi_strain_isa_c`, `materi_strain_isa_eacc`, and
`group_materi_plasti_hypo_strain_isa` are registered so the Professional syntax
input can be parsed and its declared DOFs/data can be stored. The group record is
an indexed two-double MATERI record; the syntax fixture stores `20.0 0.017`, and
no constitutive semantics are inferred from those values. The GNU incremental
driver does not read any ISA record. This is intentionally not ISA/Niemunis
constitutive support: no ISA state evolution, material dispatch, or execution
behavior has been added. See `incremental_driver_isa_parse_only.md` for the
measured Professional layout.
