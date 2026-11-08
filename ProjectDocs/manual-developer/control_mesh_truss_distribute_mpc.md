# control_mesh_truss_distribute_mpc

`control_mesh_truss_distribute_mpc` and `_air`, `_dof`,
`_element_group_truss` / `_isoparametric`, `_exact`,
`_exact_minimal_length` (`_connect`), `_geometry_truss` /
`_isoparametric` (manual Professional 6.245-6.254; F5 T4 U7).

## Files and functions

- `mpc.cc`
  - `mesh_truss_distribute_mpc()`: entry point, called from `step_start()`
    (`top.cc`) after `extrude()` and before `interface_convert()` /
    `data()` / `merge()`. A static flag makes it run once, at the first
    processed control index that has the record (also an index without
    `control_timestep`, as the Professional `other/truss10`). More than one
    active index is an error. Steps:
    1. reads the switches and the dof list (`mpc_dof_list()`, default the
       principal dofs; an unknown dof is a `db_error`);
    2. `td_selection()` reads the truss side and the isoparametric side
       (group list or geometry pairs, both on one side is the Professional
       error text); unequal counts on both sides is an error; the number
       of passes is the number of truss selectors (1 without);
    3. collects the isoparametric candidates (`td_iso`: element, name,
       nodes, `node_start_refined` coordinates, bounding box) for which
       `td_corners()` knows the element (dimension `ndim`, no interface
       group) and builds a per-pass mask (`td_element_selected()`, geometry
       through `element_in_geometry()` with `geometry_method_get()`) and a
       per-pass flag of the nodes of the masked elements;
    4. per pass and truss (`-truss`, 2 nodes, not processed yet): the
       intervals of the truss parameter inside each masked element
       (`td_clip()`), sorted by entry parameter; the split parameters (only
       with `_exact`: the entries of the intervals after the first); the
       pieces, dropped when shorter than the minimal length (with
       `_connect -yes` and a kept predecessor the next piece starts at the
       dropped piece's start: `alias[]`) or, with `_air -no`, when the piece
       centre is in no interval; the split nodes are created on use with
       `create_node()` (coordinates and `node_dof` interpolated between the
       truss nodes); every kept piece becomes a new element with
       `create_element()` (copies the element data of the old truss), the
       old truss is deleted with `delete_element()`;
    5. deletes the old truss nodes no longer used by any element;
    6. writes the ties after the highest `mpc_node_number` index with
       `td_tie_node()`: first the original truss nodes that lie in an
       element, then the nodes of the new trusses in element order (the
       Professional record order); a node of a masked element is not tied;
       an air node of an `_exact` piece whose centre lies in an element is
       tied at the clamped parameter (the exit point) of that element;
    7. `mesh_has_changed()`.
  - `td_clip()`: Cyrus-Beck clipping of the segment against the corner
    sides (2D: edges of the corner polygon; 3D: corner faces with a Newell
    normal and the face centroid), oriented outward by the element centroid.
    1D: the interval of the end nodes. A segment that only touches an
    element (interval <= `TD_EPS_T`) misses it.
  - `td_tie_node()`: `point_el()` with `MPC_EPS_ISO` and
    `mpc_weights_polish()` (machine-precision factors), then
    `mpc_tie_put()`.
  - `mpc_tie_put()` / `mpc_dof_list()`: the record writer shared with
    `mpc_element_group_index()` and `mpc_geometry_generate()` (zero-weight
    masters dropped, optional `mpc_from_*` marker, optional factor record).
- `database.cc` — the records (CONTROL class, `data_required` = the parent
  record, `_exact_minimal_length(_connect)` require `_exact`); the group
  and geometry lists are `DATA_ITEM_SIZE`, variable length (before U7
  `_element_group_truss` had a fixed length 1); marker
  `mpc_from_truss_distribute_mpc` (INTEGER, 1).
- `tochnog.h` / `tochnog-mod.h` — enums after
  `CONTROL_MESH_TRUSS_DISTRIBUTE_MPC_ELEMENT_GROUP_TRUSS`; prototype in
  `tochnog.h`.
- `scripts/build_safe.sh` — `td_check` (trusses with their groups and the
  normalized mpc records) and the error tests.

## Implementation details

- Hard-coded: `TD_EPS_T` 1.e-10 (parameter tolerance), the default minimal
  length `TD_MINIMAL_LENGTH_FACTOR` 1.e-11 times the truss length (the
  Professional keeps a 1.e-10 piece of a truss of length 8 and drops a
  1.e-12 one; the same ratio at a 1000 times larger scale), the clipping
  tolerance 1.e-10 times the element size, the bounding-box margin 1.e-6
  times its size, `MPC_EPS_ISO` (1.e-4) for the node location and 10 times
  it for the exit point of an air node.
- Supported isoparametric elements: bar2, bar3, tria3, tria6, quad4, quad9,
  tet4, tet10, hex8, hex27, prism6 (straight sides; a curved quadratic side
  is taken through its corners). Other element types are not candidates.
- The coordinates are `node_start_refined`.
- The new node keeps the node data of the first truss node (`create_node()`
  copies items present in both truss nodes); `node_dof` and
  `node_dof_start_refined` are interpolated.
- Verification: the Professional 02-08-2026 examples `other/truss10`,
  `truss11`, `truss12`, `truss13`, `truss15` (copies in test-2026) and
  about 50 probes (1D bars, 2D quad4, 3D hex8, air, minimal length,
  connect, dof, groups, geometries, errors); each test-2026 `.dbs` equals
  the Professional one after matching the nodes by coordinates, except the
  manual cases. See SEGUIMIENTO-CONVERGENCIA.md, row F5 T4 U7.

## Pending refactorings

- The candidate loop is O(trusses x elements) and the node location
  O(nodes x elements); a spatial search would help large meshes.
- Curved quadratic sides are clipped as straight ones.
- `-trussbeam` generation (manual "special option", no record defined) is
  not implemented.
- The mpc solve converges slowly for a stiff truss tied across an interface
  with the default 2 iterations (BACKLOG C128, pre-existing in the mpc
  solve, not in the generation).
- The run-once guard is a static flag: a restart from a `.dbs` that still
  holds the record distributes the already distributed trusses again (new
  numbers, duplicated ties); the generated `mpc_from_truss_distribute_mpc`
  markers could be used to skip it. Not measured on the Professional.

## Differences with the Professional version

Listed in DIVERGENCES-PROFESSIONAL.md: TRUSS-DISTRIBUTE-AIR-DEFAULT,
TRUSS-DISTRIBUTE-AIR-NO, TRUSS-DISTRIBUTE-AIR-EXIT and
TRUSS-DISTRIBUTE-ISO-SELECTION (class D: the GNU follows the manual),
TRUSS-DISTRIBUTE-NUMBERING (class C: compact numbering, no zero-factor
masters), C128-MPC-TRUSS-CONVERGENCE (class B, BACKLOG C128). `_dof` takes the control index as the Professional although the
manual heading shows none (no divergence with the Professional).
