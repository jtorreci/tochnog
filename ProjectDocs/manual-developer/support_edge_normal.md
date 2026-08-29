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
