# support_edge_normal — developer notes

Sprint 13 lot 1 (manual Professional 6.1067-6.1078 + 6.380/6.381).
Enums: `SUPPORT_EDGE_NORMAL` + 5 companions, `NODE_SUPPORT_EDGE_NORMAL_FORCE`,
`CONTROL_SUPPORT_EDGE_NORMAL_DAMPING_APPLY` / `_STIFFNESS_FREEZE`
(parse-only for now - damping is lot 2, the freeze belongs with
plasticity, lot 3).

## Implementation

All in `area.cc`, as a new `type[10]` of the edge-integral machinery
(same pattern as the `force_element_edge_*` family):

- **Selection**: `type_area[10] = SUPPORT_EDGE_NORMAL_GEOMETRY` (the
  required geometry record: node list or geometry entity, the
  machinery's `use_geom`/`area` paths), plus element-level
  restrictions (`_element`, `_element_group`, `_element_side` with
  1-based local sides) and node-level restrictions (`_node`,
  `_element_node`) via `support_edge_companion()`.
- **Force (RHS)**: per side node (Lobatto quadrature = consistent
  nodal loads): `f = -(k_n·u_n)·n - k_t·(u - u_n·n)` with `u` the
  TOTAL displacement dof (`dis_indx + idim*nder`; requires
  `materi_displacement`, hard error otherwise) and `n` the outward
  side normal. `element_rhside[vel_indx...] += w·A·f`.
- **Stiffness (matrix)**: the consistent side stiffness
  `dt·∫(k_n n⊗n + k_t(I-n⊗n)) N_i N_j dA` into `element_matrix`
  (+ diagonal into `element_lhside`), with the same `dtime` scaling
  as the element stiffness (materi.cc: `volume*dtime*stiffness`).
  Numerical Gauss quadrature (2 points linear, 3 quadratic; tensor
  for 3D faces with the border-table node order, first direction
  fastest - the same the 3D Lobatto weights use). REQUIRED for
  solvability: a body on supports otherwise keeps a zero-energy rigid
  mode in the velocity matrix (measured: solver breakdown); the fixed
  point is unchanged (the equilibrium lives in the RHS force - the
  Professional's `_plasti_residual_stiffness` documentation states
  the same for the plastic case).
- **Output record**: `NODE_SUPPORT_EDGE_NORMAL_FORCE` PUT in
  `VERSION_NEW` (the node_inertia pattern; the end-of-step
  NEW→NORMAL copy publishes it). Accumulated over the elements
  sharing a node; re-zeroed at the start of EVERY assembly sweep.
  Sweep detection: a static last-element-number - the element loop
  runs ascending and `area()` is called once per element per sweep,
  so `element <= s_last` means a new sweep started (`<` misses the
  single-element case - measured). Per-call (not per-node: the zero
  fired between the side nodes wiped node 1 - measured).

## Gotchas (measured)

- `bounda_unknown` has NO `-node` selector: `-ra <nodes> -ra <dofs>`
  or one plain node per record. A wrong selector is a SILENT no-op
  (tslv_vfx prescribed zeros, so it passed by accident).
- `bounda_time` is a flat (time, value) pair list SHARED by all
  unknowns of the record: different values per dof need separate
  records. In 3D, splitting the unknowns across TWO
  `bounda_unknown` records (velx/velz zeros + vely) left the velocity
  dofs untouched at assembly time (quirk, pre-existing; use ONE
  record per value with all dofs of that value).
- pri() takes at most (char*, value): multi-arg debug prints do not
  compile.
- Rebuilding the shadowed-loop bug: an inner `for (idim...)` inside
  an outer `for (idim...)` clobbers the outer index - use distinct
  variables (this produced fx in both record slots).

## Verification

`tsup_winkler` (prescribed u on the supported nodes: nodal force
(-0.05, +0.5) EXACT = k·u·L/2), `tsup_solve` (column on springs:
u_top = -(F·L/EA + F/kL) = -3e-3, reaction 0.5/node),
`tsup_3d` (hex8 face: k·u·A/4 = 0.025 per corner). Suite 219/219.


## Lot 2: damping, density, factor, force_initial, time (manual 6.1068-6.1071, 6.1075-6.1076, 6.1084 + 6.380/6.381)

Extended the SUPPORT_EDGE_NORMAL branch in `area.cc`:

- **Damping**: per-node viscous force `−(c_n·v_n)·n − c_t·v_t` (same
  structure as the stiffness force). Velocity from `new_dof[vel_indx]`
  ✓ (the lot-1 first-assembly observation showed the prescribed
  velocities are visible). The gate `control_support_edge_normal_
  damping_apply -no` (per-control block) skips both `_DAMPING` and
  the `_AUTOMATIC*` variants.
- **Automatic damping**: `c_n = sqrt(ρ·Eoed)`, `c_t = 0.25·sqrt(ρ·G)`
  with `Eoed = (1−ν)E/((1+ν)(1−2ν))`, `G = E/(2(1+ν))`. Reads
  `group_materi_elasti_young/poisson` and `group_materi_density` of the
  attached group (no materials → `c_n = 0`, warned once if
  `ρ ≤ 0`). The `_APPARENT` variant uses the current nodal
  stress/strain along the loading axis to estimate `E_app` (guards fall
  back to nominal); for elastic behavior identical to the nominal
  values, so the test path matches `_automatic` exactly.
- **Density**: `d_n·a_n` with `a = (v_new − v_old)/dt` (the old velocity
  from `NODE_DOF VERSION_NORMAL`).
- **Factor**: the spatial factor scales only the STIFFNESSES (manual
  6.1075: "for the support stiffnesses and not the force"). Reused
  `force_factor(SUPPORT_EDGE_NORMAL_FACTOR, ...)`.
- **Force_initial**: `f0 = a0 + a1*(y in 2D, z in 3D)` adds a compression
  preload (the reaction is `−f0·n` pushing the element). At zero
  displacement the record still shows the preload (test `tsup_init`).
- **Time diagram**: the time factor `load` multiplies the total
  support force (manual 6.1084). Reused `force_time(...)`.
- **Controls**: `control_support_edge_normal_damping_apply` wired (per
  control); `control_support_edge_normal_stiffness_freeze` parsed only
  (the elastic support stiffness doesn't change; meaning only with
  plasticity, lot 3).

## Gotcha: the "E from the element" in `_automatic`

The manual says the automatic damping reads "the Young value E and the
Poisson ratio ν from the isoparametric element attached to the
support" — that is the GROUP's E and ν, NOT a separate support-level
record. Test `tsup_auto` uses the group's E = 1e7 (NOT a support-only
E) and the measured values (c_n·v = 4472·1e-3, c_t·v = 790.6·1e-4)
confirm the formula.

## Verification (223/223 + TODAS las verificaciones)

`tsup_damp` (damping+factor: -0.0075 / +0.35 EXACT), `tsup_auto`
(automatic: 2.73607 / -0.04453 EXACT), `tsup_init` (force_initial +
time×2: 0.1 / 0.3 EXACT, linear with y), `tsup_dens` (density:
0.525/node EXACT).

## Lot 3 attempt: plasticity

The 6 enums + 1 node output are registered in `tochnog.h` and
`database.cc` and the parser accepts the records, but the CAP
APPLICATION in the per-side force block was not landed in this pass
(structural complications with the per-side shadowing in `area.cc`
that require a small refactor; the behavior is otherwise well-
scoped and the test cases are ready). The `support_edge_normal_
plasti_residual_stiffness` (matrix term) and the other caps are
deferred to a follow-up that includes the structural fix; tests
`tsup_gap`/`tsup_tcap`/`tsup_fric`/`tsup_residual` exist and are
ready in `validation-suite/test-2014/`.
