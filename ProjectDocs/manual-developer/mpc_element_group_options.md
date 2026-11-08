# mpc_element_group options

`mpc_element_group_always` / `_closest` / `_coord_geometry` / `_eps_iso` /
`_keep` and `control_mpc_element_group` (manual Professional 6.256,
6.859-6.865; F5 T4 U6).

## Files and functions

- `mpc.cc`
  - `mpc_element_group_generate()`: called from `mpc_node_apply()` every
    equilibrium iteration; returns while the stored fingerprint
    (`MPC_ELEMENT_GROUP_MESH_FINGERPRINT`) equals the current one. With
    `CONTROL_MPC_ELEMENT_GROUP -yes` at the current `ICONTROL` the
    fingerprint is the mesh hash xor a hash of the `TIME_CURRENT` bits, so
    the records are created again once per timestep. On regeneration it
    deletes the previous generated range, drops kept records whose nodes no
    longer exist, then generates in two passes: first the `_keep -yes`
    indices without a `MPC_ELEMENT_GROUP_KEPT` record (their records stay
    outside the regenerated range), then the other indices.
  - `mpc_element_group_index()`: the per-index generation (the former loop
    body): reads `_always`, `_closest`, `_coord_geometry` and `_eps_iso`;
    marks the `element_group_1` nodes and their lowest last coordinate when
    `_always -no` or `_coord_geometry -yes` needs them; searches a master
    with `point_el(..., eps_iso)` (pass 0) and, with `_closest -yes`, the
    master with the smallest `mpc_closest_weights()` distance (pass 1, with
    a bounding-box lower bound to skip far elements).
  - `mpc_closest_weights()`: `mpc_closest_tensor()` for bar / quad / hex
    (every face of the isoparametric box, each coordinate free / -1 / +1,
    solved by Gauss-Newton over the free coordinates; the Euclidean closest
    point) and `mpc_closest_simplex()` for tria3 / tria6 / tet4 / tet10
    (barycentric coordinates on the corners, the non-first ones clamped to
    [0,1], a negative first one set to 0 and the rest normalized: the
    Professional tria3 rule). Prisms and other elements are not candidates.
- `database.cc` — the five per-index records (INTEGER / DOUBLE_PRECISION,
  length 1, CONTROL, `data_required = MPC_ELEMENT_GROUP`),
  `control_mpc_element_group` (INTEGER, 1, CONTROL) and the internal
  `mpc_element_group_kept` (INTEGER, fixed length 2, `external = 0`:
  start index and count of the kept `mpc_node_number` records).
- `tochnog.h` / `tochnog-mod.h` — enums after `MPC_ELEMENT_GROUP_GEOMETRY`
  and `CONTROL_MPC_APPLY`.
- `scripts/build_safe.sh` — `mpc_check` (normalized `mpc_node_number` /
  `mpc_node_factor` records) and the fourteen tests.

## Implementation details

- Defaults: `_always -yes` (the previous GNU behaviour: master elements that
  contain the node are skipped, the search continues), `_closest -no`,
  `_coord_geometry -no`, `_eps_iso` `MPC_EPS_ISO` = 1.e-4, `_keep -no`,
  `control_mpc_element_group -no`. Without the new records the generated
  records are the same as before (differential of 35 inputs with mpc
  records: identical).
- `_eps_iso < 0` is a data error (`db_error`).
- `_coord_geometry`: the comparison is `coord[ndim-1] < lowest` (a node at
  the lowest level is tied, as in the Professional); no `element_group_1`
  node gives `lowest = DBL_MAX`, so no node is tied.
- `_keep` records live outside the regenerated range, so `mpc_node_apply()`
  treats them as user records; the condensation rule is the same for the
  generated range and user records (`in_eg_range || from_bounda!=-YES`).
- The records are read with explicit record names (the buffer audit cannot
  follow a record number held in a variable).
- Verification: Professional 02-08-2026 probes (no corpus example uses the
  records); see SEGUIMIENTO-CONVERGENCIA.md, row F5 T4 U6.

## Pending refactorings

- The factors of a node accepted outside its master within `_eps_iso` are
  extrapolated by `mpc_weights_polish()`, and `point_el()` misses such nodes
  near a corner (BACKLOG C127); clamping them as `mpc_closest_weights()`
  does would match the Professional.
- `mpc_element_group_index()` scans all elements for every node (O(nodes x
  elements)), as before; a spatial search would help large meshes.
- The hard-coded `MPC_WEIGHT_ZERO` (1.e-6) and the Gauss-Newton limits (50
  iterations, 1.e-15) of the closest point search.

## Differences with the Professional version

Listed in DIVERGENCES-PROFESSIONAL.md: MPC-ELEMENT-GROUP-ALWAYS,
MPC-ELEMENT-GROUP-COORD-GEOMETRY, MPC-ELEMENT-GROUP-EPS-ISO,
MPC-ELEMENT-GROUP-KEEP, CONTROL-MPC-ELEMENT-GROUP (class D: the GNU follows
the manual), MPC-ELEMENT-GROUP-CLOSEST (class C) and C127-MPC-OUTSIDE-WEIGHTS
(class B). The Professional also numbers its generated records from the
`mpc_element_group` index, lists zero-factor masters and prints an
`mpc_from_element_group` record (pre-existing, not part of U6).
