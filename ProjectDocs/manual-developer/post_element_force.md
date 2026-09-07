# post_element_force

Developer manual of the `post_element_force` family (Professional manual
6.927-6.937): the section resultants (normal force, shear forces and
bending moments) computed from the ELEMENT INTERNAL NODAL FORCES
(`f_elem = int B^T sigma dV`, the same free-body statics as the
`materi_stress_force` support machinery; the helpers
`msf_element_rule` + `msf_element_internal_forces_2d/3d` are shared).

## Implementación (2D/3D, commits 53ba151 + 805c445)

- Record family: `post_element_force` (the direction frame: normal +
  shear0 + shear1 + the middle point; 3D 12 values, 2D 7, 1D 2) +
  `_geometry` (nodes on a geometry item, INITIAL location) + `_group` +
  `_number` + `_normal` + `_force` + `_inertia` + `_multiply_factor`;
  result record `post_element_force_result` (5 values: normal,
  shear0, shear1, moment0, moment1).
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
