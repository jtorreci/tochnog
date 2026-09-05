# control_mesh_cut_geometry

## Files and functions

- `delete.cc` — `void mesh_cut( double time_current )` (next to `delete_geom`): the whole cut machinery. Static helpers in the same file:
  - `mc_element_rule` — the integration rule of one element (replica of the `msf_element_rule` helper of `calcul_force.cc` / `pol()` in `polynom.cc`): `-bar2` integrates with the MINIMAL 1-point rule (its stiffness is exact with 1 point, `polynom.cc:397-398`), the other tensor-product elements with the MAXIMAL rule of their group; SRI switches to Gauss.
  - `mc_element_internal_forces` — the element internal forces `f_elem = int B^T*sigma dV` at the element nodes from the `ELEMENT_DOF` integration-point stresses of the converged state (the same quantity `materi()` builds into the element right-hand side with the opposite sign; same source and kinematics as the validated section-force integration of `calcul_force.cc` LOT 5).
- `top.cc` — called in `step_start` right after `delete_geom(time_current)`; early-returns when no `CONTROL_MESH_CUT_GEOMETRY` record is active at the current control index.
- `geometry.cc` — the `CONTROL_MESH_CUT_GEOMETRY` projection type is treated exactly like `CONTROL_MESH_DELETE_GEOMETRY` in the 8 `geometry()` branches that select the "filled interior" membership for area geometries (circles, spheres, cylinders, ...).
- `database.cc` — registrations + the alias `control_mesh_cut_force` → `CONTROL_MESH_CUT_NODE_FORCE` in `db_number` (the corpus/Professional record name is `control_mesh_cut_node_force`; the manual section 6.164 writes `control_mesh_cut_force`).
- `tochnog.h` / `tochnog-mod.h` — enums `CONTROL_MESH_CUT_GEOMETRY`, `CONTROL_MESH_CUT_NODE_FORCE` appended before `LAST_DUMMY` (both headers, same order).

## Implementation details

Semantics measured on the Professional binary 25-10-2023 (`mesh_cut_1.dat` / `mesh_cut_2.dat` .dbs + stdout, plus a 2D area-cut probe):

1. **Activation** — `mesh_cut()` runs at every `step_start` whose control index has the record active. `mesh_cut_1.dat`/`mesh_cut_2.dat` put the cut records at an index WITHOUT a `control_timestep` (index 30 between the loading block 20 and the post-cut block 40): the control loop of `top()` gives such indices their own virtual `step_start`/`step_close` pair (`top.cc` else-branch of the timestep), so the cut fires exactly once, at the transition time, from the converged state of the previous block.
2. **Membership** — the nodes inside the geometry are found with `geometry( ..., NODE_START_REFINED, CONTROL_MESH_CUT_GEOMETRY, ... )`; the elements deleted are those whose nodes ALL lie inside (the same test `delete_geom` uses for `control_mesh_delete_geometry`). Verified: in `mesh_cut_1.dat` the Professional deletes exactly the bar elements fully inside the cut line (elements 6-9, nodes 7-10 pruned, node 6 kept as the boundary node); in `mesh_cut_2.dat` NO node lies within the tolerance of the upper-half line, so the Professional deletes nothing (no-op; reproduced).
3. **Nodal-force substitution** — for each element to be deleted, `mc_element_internal_forces` integrates `f_elem` from the `ELEMENT_DOF` stress of the last converged step (the state at the cut) and the code accumulates `node_force_acc -= f_elem` over the element nodes. The interior contributions of the deleted block cancel pairwise, leaving exactly the force the removed material exerted on the surviving boundary nodes (Newton's third law: the action of the removed part = `-f_elem`).
4. **Deletion and pruning** — `delete_element` for each selected element; then `mesh_has_changed( VERSION_NORMAL )` (guarded by a `VERSION_TMP` copy like `delete_geom`) whose `nod_nod()` prunes the nodes left with no element.
5. **`NODE_FORCE` records** — written on the surviving nodes with a non-zero accumulated force, per direction switch of `control_mesh_cut_node_force` (default all `-yes`). Accumulated onto any pre-existing user `node_force`.

## Verification against the Professional

- `mesh_cut_1` (1D bar2, real cut): remaining mesh nodes 1-6 at x = 0..5/9 (same as the Professional .dbs), post-cut `post_point_dof` disx 0.25 / sigxx 1.0 sustained (target ±1e-3, rc=0), substituted boundary force magnitude 1 (Professional: `node_force 6 +1`).
- `mesh_cut_2` (2D quad4, cut line that no node lies within): nothing deleted, model continues unchanged (rc=0) — identical to the Professional (its .dbs keeps all 81 elements).
- 2D area-cut probe (quadrilateral over the top rows, pull + cut): substituted forces on the boundary row = 1/9 (interior nodes) / 1/18 (edge nodes) — the consistent nodal loads of the uniform sigma_yy = 1 — EXACTLY the Professional values (sum 1.0), remaining 45 elements, post-cut sigyy at the lower point = 1.0 (rc=0).

## Hardcoded parameters / pending refactorings

- **`node_force` sign divergence (PENDING, measured)**: the Professional applies a stored `node_force +F` as a `+x` external load; `parallel_new_dof_before` (`dof.cc`) applies it as `-x` (`node_rhside -= node_force`). `mesh_cut` therefore stores the NEGATED substitute so the equilibrium matches the Professional (mesh_cut_1: `node_force 6 = -1` here vs `+1` there; both keep sigxx = +1). A dedicated work unit should align the `node_force` record sign with the Professional (A/B against the Pro binary) and drop the negation.
- The integration replicates `pol()`/`materi()` for the tensor-product family (bar/quad/hex up to cubic). Elements outside that family abort with a clear message instead of deleting without the equilibrium substitution. Axisymmetric groups are rejected.
- `mc_element_internal_forces` uses the reference coordinates (`NODE` at the cut time); models whose mesh convects (velocity-only, default `options_mesh -follow_material`) would integrate over the moved geometry — the mesh_cut corpus family is `-total_linear` (see the relax family developer manual for the velocity-only kinematics divergence).
