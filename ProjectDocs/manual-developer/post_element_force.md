# post_element_force

Developer manual of the `post_element_force` family (Professional manual
6.927-6.937): the section resultants (normal force, shear forces and
bending moments) computed from the ELEMENT INTERNAL NODAL FORCES
(`f_elem = int B^T sigma dV`, the same free-body statics as the
`materi_stress_force` support machinery; the helpers
`msf_element_rule` + `msf_element_internal_forces_2d/3d` are shared).

## Implementación (2D/3D, commits 53ba151 + 805c445)

- Record family: `post_element_force` (the direction frame: normal +
  shear0 + shear1 + the middle point; 3D 12 values, 2D 6, 1D 2; the 2D
  length was 7 before C84) +
  `_geometry` (nodes on a geometry item, INITIAL location) + `_group` +
  `_number` + `_normal` + `_force` + `_inertia` + `_multiply_factor`;
  result record `post_element_force_result` (3D 5 values: normal,
  shear0, shear1, moment0, moment1; 2D 3 values since C84: normal,
  shear0, moment0).
- `post_element_force_calculate()` (calcul_force.cc, called in the step
  close of top.cc) scans the restricted elements, accumulates
  `node_force[node] -= f_elem` (the NEGATIVE of the internal force =
  the force the element exerts on the section nodes) and sums the
  geometry-restricted nodes: normal = sum(f.n),
  shear = sum(f.s), moment = sum(f.n)*(r-mid).s. Sign calibrated
  against the Professional force10 (3D, compression -123.4 etc.).

## The 1D sections (2026-09-07, dev/segv: the force16/17, mpc8/9,
## post7 SEGV family)

### Root cause of the 5 SEGV (rc=139)

All five corpus SEGVs (force16, force17, mpc8, mpc9, post7 - all 1D
bar models with `post_element_force`) crashed in the SAME instruction:
`post_element_force_calculate()`, calcul_force.cc:2396 (the
ELEMENT_DOF integration-point stress read), on a WRITE overflow of the
stack array `sig_ip[6*MPOINT]` driven by a garbage integration-point
count. The chain:

1. `msf_element_rule()` only filled `nper[idim]` for `idim<ndim`. For
   1D problems (ndim=1) `nper[1]`/`nper[2]` stayed UNINITIALIZED stack
   garbage.
2. `npoint_ip = nper[0]*nper[1]` then multiplied the 1D point count by
   garbage -> a huge loop count -> `sig_ip[ip*6+c]` wrote far beyond
   the 162-double array (up the stack, past the guard page for mpc8 -
   its death was silent because the kernel kills the process before any
   SIGSEGV handler runs).

Verified: all five backtraces end at the same offset
(calcul_force.cc:2396); removing the MPC records does not change the
crash (they only exercise the same post_element_force code); under ASan
the garbage happened to be small/zero and the runs completed (mpc8:
"actual value 0") - the classic uninitialized-stack signature.

### Fix (calcul_force.cc)

1. `msf_element_rule()`: initialize `nper[0]=nper[1]=nper[2]=1` (the
   untouched directions of a 1D element read 1, mirroring pol(), which
   integrates the directions idim>=ndim with a single point).
2. `msf_element_rule()`: default the -bar2 integration to the MINIMAL
   1-point rule (pol() polynom.cc:551-552 integrates the -bar2 with
   MINIMAL; the other tensor-product elements with MAXIMAL). The
   section integration must replicate the rule the element really
   integrated with or the ELEMENT_DOF stresses integrate to the wrong
   force (the MAXIMAL 2-point Lobatto rule reads an all-zero second
   ELEMENT_DOF block -> half force; same gotcha as the mesh_cut
   family, delete.cc mc_element_rule).
3. `post_element_force_calculate()`: dedicated 1D branch (ndim==1).

### The 1D section semantics (measured against the Professional
### 25-10-2023 on force16/17, mpc8/9, post7 + a static probe)

The 1D result (manual 6.929/6.937: partial record `dir_normal_x
middle_x`, single normal force in `post_element_force_result`) is the
AXIAL FORCE of the element at the section (N = sigma*A, tension
positive) - NOT the 2D/3D per-node `-f_elem` sum:

- Interior section through a node shared by two elements: the -f_elem
  of both neighbours CANCEL in the node_force accumulation, while the
  Professional reports the (continuous) axial force sigma*A (post7
  middle section: -2.00, the same value as the left-end section).
- Mesh-end section: the axial force of the single adjacent element
  carries the sign of the stress (post7 fixed end: -1.81 = +f_elem of
  the last element, compression negative; force16 free tip: +5.034 =
  +sigma*A, tension positive; mpc8/9 with the group restriction: the
  axial force of the restricted element at the section node).

Per section node the implementation uses the adjacent restricted
element on the +dir_n side of the node when present (its -f_elem at
the node = its axial force) and the -dir_n side element otherwise
(+f_elem at the node = the same axial force). The axial force itself
is orientation-independent: f_elem at the element's high-coordinate
end node is +N always (the delete.cc 1D kinematics with the SIGNED
Jacobian). `post_element_force_normal -yes` keeps the 2D semantics
(only the elements on the positive side of the plane through the
middle contribute).

### Resulting corpus states (blast radius force*/mpc*/post* + suite)

- mpc8 (section x=1 through the shared node, group 2): PASS rc=0
  (0.5 within 1e-8, = sigma(el2) of the ELEMENT_DOF record).
- post7 (dynamic 1D bar, sections at x=0 / 50 / 100 +
  post_node_result at the fixed end): PASS rc=0 (all four targets
  within 1e-2).
- force16 (dynamic bar dragged by a prescribed velocity):
  RUNFAIL rc=1 - the section machinery now reports the element's own
  axial force EXACTLY (2.93667 = the ELEMENT_DOF sigma of bar 1), but
  the GNU solved sigma field differs from the Professional's
  (sigma1 = v2-v1 = 2.93667 vs 5.03429; v2 = 3.9367 vs 6.0343). Same
  family as the documented dynamic* scheme gap (GNU staggered
  u-sigma temporal map vs the Pro displacement-primary one; see the
  dynamic* closure) - not a section bug.
- force17 (no geometry + `post_element_force_force -yes` +
  `post_element_force_inertia -yes`): RUNFAIL rc=1 - the
  force/inertia switches are PARSE-ONLY in the GNU fork (registered in
  database.cc, never consumed), and the no-geometry elastic sum
  (2*sigma(el1) = -2.478) cannot reproduce the Pro's
  elastic+external+inertia value 9.93172. Documented RUNFAIL.
- mpc9 (explicit tie between two free elastic nodes of OVERLAPPING
  bars, x=0.8 tied to x=1.0): RUNFAIL rc=1 - the GNU explicit-tie
  semantic (value-constrained slave, row excluded WITHOUT force
  redistribution, mpc.cc:25-37 - measured against the Pro on mpc1)
  leaves the master's equation free of the slave's elastic share:
  v2=v3=0, sigma1=0, sigma2=0.833 vs the Pro v2=v3=0.4545,
  sigma1=sigma2=0.4545. The energy-consistent slave elimination (fix
  G) is deliberately activated ONLY for the mpc_linear_quadratic
  generated ties (mpc3/mpc4 closure); extending it to explicit user
  ties is an open semantic question with blast radius over
  mpc1/2/7/8 (byte-identical targets) - out of the SEGV scope.

### Files

- calcul_force.cc: msf_element_rule (nper init + bar2 MINIMAL) +
  the ndim==1 branch of post_element_force_calculate.

## C84: element types, the 2D record and axisymmetric groups (F4 T6f, 2026-10-17)

### Problem

The 2D/3D branch of `post_element_force_calculate()` picked `npol`
from the element name (2 for quad4/tria3/hex8, 3 otherwise), built a
tensor-product rule with `msf_element_rule()` and integrated with
`msf_element_internal_forces_2d/3d()` for every element type. A tria,
tet or prism element got quad/hex shape functions, a 2x2(x2) or
3x3(x3) point count and read `ELEMENT_DOF` slots past its real points
(the `tpef_force10` cantilever split into tet4: normal 3143.7, shear1
-539.6, moment1 -16896 instead of -123.4 / 1000 / -100000). In 2D the
record required 7 values (manual 6.927: 6), and the stress copy wrote
`sxx sxy sxz` into a 6-stride array that the 2D helper reads as `sxx
syy sxy` with stride 3: uninitialized values, NaN for every 2D result,
quad4 included. The 2D result had 5 values; the Professional writes 3.

### Fix (`calcul_force.cc`, `input.cc`, `database.cc`)

- `input.cc` (data part, before the store): `post_element_force` must
  have 12 / 6 / 2 values for ndim 3 / 2 / 1, else `Error in data part.`
  + `Number of data values expected : N` + the count found, exit 1
  (the Professional message, which also stops in the data part).
- `post_element_force_calculate()`: the 2D frame reads 6 values; an
  invalid length left by another path is an error (was a skip warning).
- Element routing in the 2D/3D element scan:
  - tensor path, unchanged: quad4/quad9 (plane) and hex8/hex27 keep
    `msf_element_rule()` + the quad/hex helpers. The 2D copy now fills
    `sig_ip[ip*3+{0,1,2}]` with `stress_indx(0,0)`, `(1,1)`, `(0,1)`
    (the helper layout); the 3D copy is unchanged (hex8/hex27 results
    byte-identical: `tpef_force10` and hex8/hex27 cantilever probes,
    masked `.dbs` diff 0).
  - pol path, new `pef_element_internal_forces_pol()`: tria3, tria6,
    quad16, tet4, tet10, prism6, prism15, and every 2D element of an
    axisymmetric group. It calls `pol()` (polynom.cc) with the node
    coordinates (ndim stride) for the real point count, the physical
    derivatives `new_d[ip*ndim*nnol + j*nnol + inol]`, `h` and the
    point volumes, and integrates `f[inol,i] += vol * sigma_ij *
    dN_inol/dx_j` with the `ELEMENT_DOF` stresses of those points.
    Axisymmetric: the pol() volume keeps 2 pi r and the radial force
    gets the hoop term `vol * sigma_zz * N/r` (the element's own ring
    forces). `check_element_shape_factor` is set to 0 around the call
    (the distortion warning was already printed by the solve); pol()
    leaves its thread-local recovery/HRZ state as for any element visit.
  - anything else (bar/truss/beam/spring names, `hex64`, an element in
    a `group_interface` group): `Error: post_element_force does not
    support element : N` + the supported list, exit 1.
- Both sides (manual 6.933, "elements on both sides will be used if
  present"): the element centroid side `(c - mid) . dir_n` picks the
  accumulator, `node_force` (>= 0) or `node_force_minus` (< 0); with
  `post_element_force_normal -yes` the negative side is skipped as
  before. A side is present when one of its elements touches a section
  node (`pef_node_ok`, the geometry restriction computed once per
  record). Result: positive side only = its resultant (unchanged,
  hex8/hex27 boundary sections byte-identical); both sides = (R+ -
  R-)/2; negative side only = -R-. Before C84 both sides were summed
  into one array, which cancels at an interior section to the
  external load at the section nodes (`force14`: shear 33.3 instead of
  the Professional -500).
- 2D result: `res[2] = moment0`, length 3; `database.cc`
  `fixed_length[POST_ELEMENT_FORCE_RESULT] = 0` (3D keeps 5 values,
  the 1D branch is unchanged).

### Verification

Professional 02-08-2026 on the same files (tolerance 4e-6 of
max(|value|,1)), `test-2026`:

- `post_element_force_{tet4,tet10,prism6,prism15}_1`: -123.4 / 0 /
  1000 / 0 / -100000 (GNU before: tet4 3143.7, tet10 -1651.5, prism6
  4918.6, prism15 1871.5 normal).
- `post_element_force_2d_{quad4,tria3,quad9,tria6}_1`: 100 / 10 / -40
  (before: 6 values skipped, 7 values NaN / 1e168).
- `post_element_force_2d_axisym_{quad4,tria3}_1`: ring strip r 1..2,
  normal 100 (= the applied ring force), shear0 -0.2696 / -1.6269,
  moment0 6.2651 / 5.9225; quad9 and tria6 probes equal too (before:
  NaN; the plane volume alone gives 10.6).
- `post_element_force_2d_length_error1`: 7 values -> the data-part
  error (expected-error test).
- `post_element_force_2d_interior_{quad4,tria6}_1`: section at x=2 of
  the strip, both sides: 100 / 10 / -20 (plain sum: 0).
  `post_element_force_tet4_reversed1`: normal 0 0 -1, negative side
  only: -123.4 / 0 / -1000 / 0 / -100000. Corpus `force14`/`force15`
  (quad9 beam with an edge load, interior section, plane and
  axisymmetric) now equal the Professional: -500 / 12500 and
  4.2573 / -235619.449 / 2224130.70 (before C84: no result, the 6-value
  record was skipped).

RED on `54b8170` (all 14), GREEN after the fix.

### 1D result length (F4 T6f LOW queue Q1, 2026-10-17)

The Professional 02-08-2026 writes ONE value in 1D (`force16.dbs`:
`post_element_force_result 10 5.034285714286`; a `target_item` on
value number 1 stops with `Error detected for data item :
target_item`). Both 1D exits of `post_element_force_calculate()` (the
`_force`/`_inertia` lumped branch and the elastic section device) now
PUT `length = 1` (`res[0] *= fac`); before they wrote the normal force
and four zeros. `database.cc` keeps `data_length 5`, `fixed_length 0`
(the maximum, 3D). No code reads values 1-4 of a 1D record. Tests:
`test-2026/post_element_force_1d_length1` (force16, normal
5.034285714286 = Professional) and `_1d_length_error1` (target on value
1 -> expected error; it passed with value 0 before). RED on `b938b0f`
(`_length_error1` rc 0), GREEN. Differential force16/17, mpc8/9,
post7, tpef_force10, force13/14 vs `b938b0f`: only the 1D record loses
its four zeros.

### `_inertia` previous-velocity snapshot (BACKLOG C25 item 1, 2026-10-17)

`pef_vel_prev` (file-static in `post_element_force_calculate()`) is
zeroed after every reallocation, not only on the first call: a node
inactive at the refresh and active later read uninitialized memory as
its previous velocity. No test: the read needs a node activated after
the snapshot grew, and uninitialized memory is not deterministic; the
change touches only those slots (active nodes are refilled from
`NODE_DOF` right after), build_safe and the omen A/B unchanged.

### Pending

- The Professional behaviour for trusses/beams/springs or interface
  groups inside the selection was not measured (GNU refuses them;
  DIVERGENCES C84-REFUSE).
- Hardcoded: the supported name lists in `post_element_force_calculate`
  (ndim 2: quad4/quad9 tensor, tria3/tria6/quad16 pol; ndim 3:
  hex8/hex27 tensor, tet4/tet10/prism6/prism15 pol).

