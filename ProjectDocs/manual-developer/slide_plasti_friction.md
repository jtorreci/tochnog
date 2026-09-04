# slide_plasti_friction — developer notes

Slide family (slide_geometry law) of the Professional (manual 6.1041-
6.1047 + 6.370-6.372 + 6.893). The previous batch registered the records
(`database.cc` slide block: SLIDE_PLASTI_FRICTION/TENSION/STIFFNESS/
RESIDUAL_STIFFNESS + CONTROL_SLIDE_PLASTI_APPLY/STIFFNESS_APPLY; enums in
`tochnog.h` + `tochnog-mod.h` in sync) and diagnosed that the legacy
`slide()` penalty (velocity) cannot reproduce the Professional inputs
(slide1/4: block falls, top reaction ~ 0).

## Implementation (this lot)

All in `slide.cc`, called from `top.cc` in the equilibrium loop after the
element assembly:

- `slide()` iterates the slide geometry indices; per member node (start
  coordinates on the geometry via `geometry(NODE_START_REFINED,
  PROJECT_EXACT)` OR the explicit `node_slide` record) it dispatches to:
- `slide_penalty_law()` — the LEGACY velocity penalty + friction path,
  kept bit-for-bit for slide geometries WITHOUT a `slide_stiffness`
  record (classic inputs, e.g. examp14). Note: the original code read
  `SLIDE_AXISYMMETRIC` with index 0; preserved.
- `slide_spring()` — the PROFESSIONAL elastic-plastic law on the TOTAL
  displacements (per node, stateless):
  - displacement source: `dis` dofs (`materi_displacement`) or `veli`
    dofs (`materi_velocity_integrated`); no displacement dofs -> no
    spring (nothing to measure against the fixed plane);
  - `un = u.n`, `S = -kn*un` clamped at `-sig_t` (two-sided spring;
    `slide_plasti_tension` caps the pull - an ABSENT record means NO
    tension limit, NOT a compression-only gap: a gap makes the
    redistribution transient of the dragged block lift the loaded node
    permanently - measured on slide4);
  - tangential elastic predictor `kt*|u_t|` capped at `c + mu*S`
    (record order phi c -> `mu = tan(values[0])`, `c = values[1]`!);
  - RHS: force ON the node `S*n + F_t` (contact.cc convention,
    `node_rhside[vel_indx/nder+idim] += ...`);
  - matrix: per-direction diagonal `dtime*|dir_i|` weights (contact.cc
    pattern). STICK: `kt`; SLIP (capped): only the residual fraction
    `res_t*kt` (default 1e-2 = the Professional default when the record
    is absent) — the ELASTIC kt in the matrix during slip biases the
    2-iteration velocity steps (measured: slide1's x reaction -0.00466
    instead of the cap -0.005773);
  - axisymmetric: `ax = 2*pi*r` multiplies forces and matrix terms;
    triggered by the `slide_axisymmetric` record (index 0, legacy) or by
    any problem group with `group_axisymmetric -yes` (corpus slide2/3).
  - Output records per node (VERSION_NEW PUT, unconditional for member
    nodes): `node_slide_direction` (n + t = -u_t/|u_t|, the friction
    direction ON the node), `node_slide_f` (`|F_t| - (c + mu*S)` when a
    friction record exists) and `node_slide_force`
    (`-S*ax`, `-|F_t|`) - the force the material applies ON the slide
    geometry, expressed in the local (n, t) frame.

### New records in database.cc

- Outputs: NODE_SLIDE_FORCE (ndim), NODE_SLIDE_F (1), NODE_SLIDE_DIRECTION
  (6) - node class, version_all=1 (enums inserted after NODE_SLIDE in
  both headers, in sync). Registered as parse-able names so the
  `target_item -node_slide_force node slot` tests of slide2 resolve.
- PRINT_DEBUG (`print_debug -yes/-no`, no_index PRINT class, parse-only):
  the corpus slide3.dat uses it.

### input.cc: `axisymmetric -yes` without an index

The Professional allows the no-index global spelling `axisymmetric
switch` (manual 6.17); the GNU data model is per-group
(`group_axisymmetric index switch`). The data-part parser now accepts the
bare spelling: when a GROUP_AXISYMMETRIC record's index token is a
switch, the index defaults to 0 and the token is re-injected as the first
data value via a one-token static pushback (`saved_first_value`).
Rescued axisym1 too.

## Law verification (Professional 25-10-2023 .dbs)

slide1 (block 1x1, E=1, top pushed down 1e-2 and dragged +x at 1e-2,
kn = kt = 1e3, phi = 0.523598):
- final spring forces per node = kn*penetration EXACTLY (per node, no
  L/2 factor: F1+F2 = 9.995e-3 = the top reaction), friction = mu*Fn at
  the cap (node_slide_f = 0);
- the moment of the drag couple redistributes the vertical load:
  S1 = 2.11e-3 (x=0) / S2 = 7.88e-3 (x=1) - GNU reproduces both to
  4-5 digits and passes the targets (+0.01 / -0.5773e-2, rc=0).

## Blockers (documented, with diagnosis)

- slide2 (axisymmetric, target tolerance 1e-9): the GNU axisym element
  splits the nodal forces of the compressed annulus EQUALLY (2*pi*r at
  the edge midpoint) while the Professional (with the lobatto
  integration of the test) splits them proportionally to the NODE radius
  (node_slide_force 2.094e-3 : 4.188e-3). PRE-EXISTING element
  discrepancy (reproduced on a fixed-bottom axisym column WITHOUT slide
  records: GNU 4.712/4.712 vs Pro 3.142/6.283 with lobatto), so the
  slide springs settle at u1 != u2 and sigma_yy = -6.447e-4 instead of
  the target -6.666666e-4. Not a slide-law issue.
- slide4 (drag at vx = 1, 100x slide1's rate): the per-step equilibrium
  of the capped slider oscillates at the drag start (fixed iteration
  budgets from 2 to 40 all fail to converge; the node at x=0 drifts into
  lift-off). slide1 (slow drag) converges exactly with the same law, so
  the blocker is the iteration dynamics of the stiff slider with the
  fast drag: needs the incremental (stateful) elastic-plastic slip with
  plastic-slip history per node + consistent tangent instead of the
  stateless total-displacement return mapping.

## Gotchas (measured)

- Record order of slide_plasti_friction is `phi c` (swapping them gives
  mu = tan(c) = 0 and a never-capping spring - slide1's x reaction then
  grows elastically past the cap).
- A compression-only normal spring (gap at un > 0) is a trap: once the
  transient of the dragged block unloads a node past zero it can never
  re-contact and the whole vertical load jumps to the other spring
  (fold/lift-off measured at drag step 2). The two-sided spring with the
  sig_t cap avoids it.
- The matrix tangential stiffness during plastic slip must NOT be the
  elastic kt: with only 2 equilibrium iterations per step it biases the
  converged state (slide1: -0.00466 vs the cap -0.005773). The residual
  fraction (matrix-only) is the safe choice.
- PUTs of the per-node output records must be UNCONDITIONAL for member
  nodes: the db_active_index() guard never fires for records that are
  never created (the slide2 node_slide_force target then dies at
  exit_tn).
