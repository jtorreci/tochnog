# mpc_element_group family

## Where implemented

- `mpc.cc` (`mpc_element_group_generate`, called from
  `mpc_node_apply`): generates `mpc_node_number`/`mpc_node_factor`
  records for the nodes of `group_0` elements located inside a
  `group_1` element.
- `mpc.cc` (`mpc_node_apply`): global/per-timestep gates `mpc_apply`
  (no index) and `control_mpc_apply` (control index) switch the whole
  mpc machinery off when `-no`.
- `database.cc`: records `mpc_element_group`, `mpc_element_group_dof`,
  `mpc_element_group_geometry`, `mpc_apply`, `control_mpc_apply`
  (INTEGER, variable/fixed per family) and the internal bookkeeping
  record `mpc_element_group_mesh_fingerprint` (external=0, no_index=1,
  layout [mesh fingerprint, start_index, count]).
- `tochnog.h` / `tochnog-mod.h`: enums in sync (same order).

## Implementation details

Generation (per `mpc_element_group` record):

1. Master groups = stored values after `group_0`, or every other group
   when `-all` is present.
2. For each element of group `group_0`, for each of its nodes (once per
   record — visited array): skip when the node is not on the
   `mpc_element_group_geometry` (geometry()/PROJECT_EXACT at
   NODE_START_REFINED) or when no restricting geometry is given.
3. The containing master element is found with `point_el()` (first
   match, iso tolerance `MPC_EPS_ISO` 1.e-4). Master elements that
   already contain the node in their node list are skipped (the meshes
   are connected there; a tie would reduce to the node itself and
   freeze its equation — the `mpc_element_group_always -no` semantics,
   manual 6.861).
4. The shape-function weights of the master element at the node
   position become the master factors (`MPC_WEIGHT_ZERO` 1.e-6
   threshold, same as the mpc_linear_quadratic generator). One
   `mpc_node_number`/`mpc_node_factor` pair is written per selected dof
   (`mpc_element_group_dof` labels resolved through `dof_label`;
   default: all principal dofs).

Regeneration: the mesh fingerprint (hash of nodes/elements) is stored
in `MPC_ELEMENT_GROUP_MESH_FINGERPRINT`; when the mesh changes the old
generated records (stored index range) are deleted and re-created after
the highest active `mpc_node_number` index (coexists with
mpc_linear_quadratic records).

Gates: `mpc_apply` (global, default -yes) and `control_mpc_apply`
(current control index) are read at the top of `mpc_node_apply`; with
`-no` the whole machinery (generation AND application) is skipped.

## Verification

- `ground19_water_under_dam`: group-2 wall nodes tied 1:1 to the
  group-3 nodes along the cut-off wall; the flow target is met (the mpc
  does not change the pore-pressure field, verified A/B against the
  Professional: identical flux with and without the mpc block).
- `mpc7` (patch test: two small quad4 tied into one large quad4 under
  unit stress): `sigxx = 1.0` within `1.e-8` — validates the
  shape-function interpolation of the generated ties (mid-edge nodes
  get 0.5/0.5 factors).
- `mpc_apply`/`control_mpc_apply` parse and gate correctly (default
  -yes keeps the ties active).

## Pending

- `mpc_element_group_always`, `_closest`, `_coord_geometry`, `_eps_iso`,
  `_keep` and `control_mpc_element_group`: not implemented (the GNU
  always uses the `-no`-like member-skip semantics and the default
  eps_iso).
- `mpc_geometry` family (6.868-6.872) and `mpc_linear_quadratic` as
  before: the geometry family is registered but not consumed.
