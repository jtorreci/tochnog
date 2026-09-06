# mpc_linear_quadratic — developer notes

Manual Professional 6.873. Implemented in `mpc.cc`:
`mpc_linear_quadratic_generate()` (static), called from
`mpc_node_apply()` every equilibrium iteration; the generation itself
runs only when the mesh fingerprint changed.

## Algorithm

1. Read the switch `MPC_LINEAR_QUADRATIC` (no_index, length 1). Off ->
   nothing.
2. Compute a mesh fingerprint (cheap hash over the element topology:
   max_node, max_element, per-element node ids). If it equals the stored
   fingerprint (internal record `MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT`,
   registered with `external = 0` so it never appears in the .dbs) ->
   nothing.
3. Delete the previously generated records (the fingerprint record
   stores `[fingerprint, start_index, count]`).
4. Mark the nodes that belong to at least one LINEAR element (bar2,
   tria3, quad4, tet4, hex8). Collect the linear elements as master
   candidates.
5. For every node of every QUADRATIC element (bar3, tria6, quad9, tet10,
   hex27) that is NOT in any linear element's node list: find the linear
   element that contains it via `point_el()` (the same inverse-isoparametric
   machinery used by `post_point`), with tolerance 1.e-4 (same default as
   `mpc_element_group_eps_iso`, 6.865). `point_el` returns the shape
   function weights of the linear element at the node position.
6. Generate one `MPC_NODE_NUMBER` + `MPC_NODE_FACTOR` pair per principal
   dof (`dof_principal >= 0`), masters = the linear element nodes with
   |weight| > 1.e-12 (zero weights filtered, matching the Professional's
   records), factors = the weights. The record layout matches the
   Professional byte for byte (verified: mpc3/mpc4/mpc5 .dbs).

## Verification (Professional 25-10-2023, .dbs)

- mpc3 (2D quad9+quad4, `control_mesh_refine_globally`): generated
  records tie the quad9 mid-edge nodes on the interface to the quad4
  edge endpoints with factor 0.5 — identical master lists and factors.
- mpc4 (3D hex27+hex8, refinement): mid-edge slaves -> 2 masters 0.5,
  mid-face slaves -> 4 masters 0.25, 48 records = 16 slaves x 3 dofs.
- mpc5 (no refinement): the single tie (node 10 -> 9, 11 with 0.5)
  matches.

## Mesh-change regeneration

The fingerprint is recomputed on every `mpc_node_apply()` call (a
cheap linear scan of the element records) and the ties are regenerated
when the mesh changed (refinement in `step_close`, deletion in
`step_start`, split, failure-deletion, ...). The generated index range
starts after the highest ACTIVE user `mpc_node_number` index
(`db_highest_index`), so user records are never clobbered.

## Hardcoded parameters

- `MPC_EPS_ISO 1.e-4` — iso-coordinate tolerance (matches
  `mpc_element_group_eps_iso` default).
- `MPC_WEIGHT_ZERO 1.e-12` — weight filter for the generated masters.

## Gotchas

- The DB is DENSE (`data_ptr = data_length * index`): generated records
  must stay at small indices; huge sentinel indices are impossible.
- The generated records live in VERSION_NORMAL and are reverted on a
  `db_version_copy(VERSION_START, VERSION_NORMAL)` rewind; the stored
  fingerprint is also reverted, so the next `mpc_node_apply()` call
  regenerates them from the reverted mesh — consistent.
- The `makefile` object list and per-file compile rules must both list
  `mpc.o`; the rules use CRLF line endings, keep them.

## Pending

- `mpc3`/`mpc4` stay RUNFAIL with the tie GENERATION verified but the
  SOLVE semantics open. Diagnosis 2026-09-07 (DIAG-SOLVE-MIXTO.md §15):
  the GNU tie consumption (`mpc_node_apply`, bounda pre-solve +
  top.cc re-sync) is the VALUE-CONSTRAINED formulation — the slave dof
  is bounded to sum(factor*master) of the PREVIOUS iteration and its
  row is dropped WITHOUT the energy-consistent redistribution of its
  internal force to the masters. On the non-conforming quadratic
  interface the tied mid-edge node carries 2/3 of the edge traction
  (corners L/6): dropping the row loses it, the interface equilibrium
  is violated and the staggered fixed point is 0.454 (mpc3 refined)
  instead of the homogeneous 1/3. The earlier hypothesis (the mixed
  σ-dofs of the tied nodes stay free) is REFUTED by measurement: with
  the slaves pinned to their EXACT homogeneous velocities the free
  field is still polluted. The Professional is exact at ONE iteration
  (implicit). The complete fix = the implicit energy-consistent slave
  elimination inside the staggered solve (prototype implemented and
  measured in 2026-09-07: fixes the fixed point to 1/3 but does not
  converge within the 2-iteration corpus default and perturbs the
  mpc7 ±1e-8 target by 2.4e-8 → reverted; family §8 C of
  DIAG-SOLVE-MIXTO.md).
- `mpc5` stays RUNFAIL: it also needs the
  `control_mesh_delete_geometry_factor` family (element half-deletion +
  the stress RESET semantics: the GNU's staggered u-sigma update
  `sigma_new = rhs/lhs` overwrites the reset value with 0 where the
  Professional keeps it — the same root cause as delete3, see
  DIAG-SOLVE-MIXTO.md).
