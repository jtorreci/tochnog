# bounda_dof

## Files and functions

- `bounda.cc` — within `bounda()`: `bounda_dof` is implemented as a functional
  alias of `bounda_unknown`.
  - Local `max_bounda_dof=0` (line 29), added to the `max_bounda` calculation:
    `db_max_index( BOUNDA_DOF, max_bounda_dof, ... )` (line 64) and
    `if ( max_bounda_dof>max_bounda ) max_bounda = max_bounda_dof;` (line 83).
  - The iboun loop treats both records as "unknown" (Dirichlet):
    `unknown = db_active_index( BOUNDA_UNKNOWN, iboun, ... ) ||
               db_active_index( BOUNDA_DOF, iboun, ... );` (lines 86–87).
  - Value read (lines 183–188): uses `BOUNDA_DOF` when it is active, otherwise
    falls back to `BOUNDA_UNKNOWN`.
- `database.cc` — keyword registration (lines 188–192): name `bounda_dof`,
  `type = INTEGER`, `data_length = MBOUNDA`, `fixed_length = 0`,
  `data_class = BOUNDA`.
- `tochnog.h` — enum `BOUNDA_DOF` (line 142), before `BOUNDA_UNKNOWN` (line 154).
- `tochnog-mod.h` — mirror enum `BOUNDA_DOF` (line 135), must stay in sync.

## Implementation details

- `bounda_dof` is a pure alias: the record is read into the same `val` array
  and processed by the identical code path as `bounda_unknown`, so it supports
  the same selection modes — node ranges (`-range`), `-all`, `-node_set` and
  geometries (`-geometry_line`, `-geometry_point`, ...) — through the shared
  `val[0]` dispatch in the node loop (lines 288–305).
- The `unknown` flag drives both the value read and the dof-label lookup; only
  one of `bounda_dof`/`bounda_unknown` or `bounda_force` may be active per
  boundary index.
- The dofs prescribed are the primary dofs (e.g. `-velx`, `-temp`, `-pres`)
  resolved through `dof_label`/`array_member`, identical to `bounda_unknown`.
- Values over time come from the matching `bounda_time` record, as for
  `bounda_unknown`.

## External dependencies

None. Reuses the existing `bounda_unknown` code path and the `geometry()`
routine for geometry-based selection.

## Hardcoded parameters / pending refactorings

- The dof lookup error `db_error( BOUNDA_UNKNOWN, iboun );` fires even when
  the input used `bounda_dof` (messages may confuse users of the
  professional keyword).
- The alias duplicates the `bounda_unknown` registration pattern in
  `database.cc`; the two records could share a single registration helper to
  avoid drift between `BOUNDA_DOF` and `BOUNDA_UNKNOWN`.
- `data_required` is unset for `BOUNDA_DOF`, so `bounda_dof` works standalone
  (like `bounda_unknown`).

## bounda_dof -topres (total pore pressure prescription)

- `db_number("topres")` resolves to the `GROUNDFLOW_PRESSURE` keyword enum
  (NOT to the dynamic dof label that `-pres` resolves to), so a `-topres`
  token stored in the BOUNDA_DOF record does not match `dof_label[]`.
- `bounda.cc` `bounda()`: the per-dof resolution loop special-cases
  `val[iu] == -GROUNDFLOW_PRESSURE` when `groundflow_pressure` is active:
  it maps the token onto the pres dof (`iuknwn = pres_indx`) and sets the
  per-record flag `topres_bounda` (reset per iboun next to `rotate`).
- The load application (`else` branch of the value setting) converts the
  prescribed total pressure `load` into the head value per node by INVERTING
  `groundflow_phreatic_coord()`:
  - `found` (phreatic level / static height covers the node):
    `new_node_dof = load - static_pressure - addtopressure`;
  - otherwise (no level): `new_node_dof = load + dens*g*z - addtopressure`
    (the no-level total is `pres_dof - rho*g*z`).
  The node coordinate is the current one (NODE + materi displacement).
- GOTCHA: the inversion uses the SAME static the machinery adds, so it is
  exact for whatever covers the node — including a
  `post_calcul_static_pressure_height` region (ground13: bottom/top rows lie
  inside the region height_ref=123, so the prescribed -20/-10 become the
  uniform head +1210, and `-to_pres` = -20/-10 / `-dy_pres` = 1210 match the
  Professional .dbs digit by digit).
- Explicit bounds always win over automatic defaults: `bounda()` zeroes
  `NODE_BOUNDED`, then `groundflow_phreatic_apply()` applies the phreatic
  conditions and finally the bounda_dof records overwrite their nodes.

## -veln: zero velocity normal to a wall (manual Professional 6.22)

`bounda_dof <index> <selector> -veln` prescribes that the selected nodes do
not move in the direction NORMAL to a plane (frictionless wall). The
normal comes from the geometry entity when the selector is a geometry
(`geometry_line`, `geometry_set`, ...) or from `bounda_normal` of the same
index for node ranges (the manual requires it there). The condition is
imposed with GENERATED mpc records — the Professional .dbs stores
`mpc_node_number`/`mpc_node_factor` + `mpc_from_bounda <k> -yes`.

Semantics measured against the Professional binary 25-10-2023 (validation_8
+ dedicated oblique/vertical wall probes at slopes -2..+2):

- one record per boundary node, ASCENDING node order (the record order of
  the Professional .dbs);
- the node normal is the FIRST-match normal of geometry() over the entities
  of a set (measured: corner node (0,1) of validation_8 lies on the
  horizontal upper_left edge (xi=1) AND on the oblique edge (xi=0); the
  entity order of geometry_set 1 decides: the horizontal edge comes first,
  so the corner is treated horizontal, while the mirror corner (1,0.9) is
  matched by the oblique edge first and treated oblique);
- slave dof = the FIRST velocity axis (x, y, z) with |n_axis| > 1e-12,
  masters = the remaining axes with |n_axis| > 1e-12, factor
  = -n_master/n_slave. Zero-component masters are OMITTED, so axis-aligned
  walls collapse to a masterless record `mpc_node_number k <node> -vely`
  (horizontal) or `-velx` (vertical) that bounds the normal dof to zero
  (mpc_node_apply accepts masterless records: value 0).
- measured records: oblique wall slope -0.1 →
  `mpc_node_number k <node> -velx <node> -vely` factor -10 (the constraint
  velx = -10*vely <=> n.v = 0 with n = (0.1,1)/sqrt(1.01)); slope +0.1 →
  factor +10, +0.5 → +2, +2 → +0.5, -1.5 → -0.667, -2 → -0.5; vertical →
  `-velx` alone; horizontal → `-vely` alone.
- a generated tie whose slave dof an EXPLICIT bounda record of the pass
  already bounded is INERT (the direct prescription wins; measured at the
  inflow corner (0,2) of the wall: the direct velx=1 bound is kept while the
  tie record is present but never applied);
- the free ties are ELIMINATED inside the solve (registered in the
  mpc_tie_* map with the mpc_linear_quadratic machinery): the slave row is
  redistributed to the master, so the wall constraint enters the same solve
  that produces the masters (a value-constrained-only slave lags one
  iteration and DIVERGES on the staggered u-sigma iteration — measured).

### Implementation (bounda.cc / mpc.cc / database.cc)

- `VELN` keyword + `MPC_FROM_BOUNDA` marker record + INTERNAL indexed
  `BOUNDA_VELN_MESH_FINGERPRINT` ([mesh fingerprint, start index, count] per
  bounda record) — enum append before LAST_DUMMY (tochnog.h +
  tochnog-mod.h), registrations in database.cc (VELN = name only; the
  keyword is validated in check.cc: materi_velocity required).
- `bounda.cc::bounda_veln_mpc()` (static, before bounda()): called per
  bounda record whose dof list (iu_start..iu_end) contains -VELN, inside
  the found-block before the node loop. Fingerprint-guarded regeneration;
  generated records start after the highest active mpc record. Node
  iteration mirrors the selector logic of bounda() (geometry:
  nodes_in_geometry + geometry() normal; ranges: bounda_normal_vec of the
  record). Per-node normal & slave/master/factor logic as above; writes
  MPC_NODE_NUMBER (+ MPC_NODE_FACTOR when nmaster>0) + MPC_FROM_BOUNDA
  (-yes) at each generated index.
- bounda() dof loop: `if (unknown && val[iu]==-VELN) continue;` — the token
  prescribes no dof by itself (the -veln record's own geometry nodes get
  nothing from the regular path; the OTHER records of the model apply as
  usual, which is what makes the direct prescription win at shared nodes).
- mpc.cc mpc_node_apply():
  - masterless records (length 2) accepted: the slave is bounded to 0 (no
    db_error; the -veln axis-aligned walls generate those);
  - records marked MPC_FROM_BOUNDA: value-apply SKIPPED when the slave dof
    is already bounded (direct prescription wins); otherwise applied and
    registered for the energy-consistent elimination (the tie map) when
    nmaster>0 and the slave is free. The applied slaves are collected in a
    static list and the post-solve re-sync (top.cc mpc_node_apply) clears
    their bounded flags first, so the free ties refresh from the freshly
    solved masters while the direct-bound ones keep their value.
- top.cc post-solve re-sync unchanged (mpc_node_apply call already there).

### Verification

- ps01/ps01_1 probes (oblique wall slope +0.1 over a quad): generated
  records identical to the Professional .dbs; steady velocity along the
  wall satisfies velx = +10*vely EXACTLY (0.915194399803/0.0915194399803
  after 20 steps vs Professional 0.9087/0.09087, ~1% transient agreement).
- validation_8 (default Bi-CG solver): rc=0, completes the 2000 steps; the
  50 generated mpc records match the Professional .dbs record-for-record;
  final wall velocities within ~1% of the Professional (inflow velx=1 kept
  at the corners; the oblique wall nodes slide with the exact tie ratio).
  NOTE: from the step where the plastic zone builds up the Bi-CG breaks
  down and the solve falls back to the direct LU retry (log spam, ~2:45
  total); with `options_solver -matrix_lapack` the run is clean (~2:00).
  This is the pre-existing staggered u-sigma/iterative-solver family
  behavior (validation_1, mpc3/4 documentation), not the -veln machinery.
