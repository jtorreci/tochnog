# post_solve_material_state_pass

Status: **IMPLEMENTED 2026-10-10** (F4 T5, BACKLOG C35, audit K02;
project-owner decision 2026-10-10). This page replaces
`nodal_stress_post_solve_sync` for every model the gate does not exclude.

## Problem

`top.cc` `top()` iterates `element_loop(); solve(); parallel_new_dof_diagonal`.
The element loop of iteration k evaluates the material on the PRE-solve
iterate (velocity of iteration k-1). It writes the integration-point state
(ELEMENT_DOF VERSION_NEW: stress, strains, plastic strain, kappa, history
variables) and assembles the non-principal nodal right-hand sides that
`parallel_new_dof_diagonal` applies after the solve. After the last
iteration the stored state and the non-principal nodal dofs therefore lag
the solved velocity by one iteration. They are the step output and the
next step's VERSION_NORMAL history. The Professional evaluates the material
once more on the final iterate: `spring4` with 1 iteration gives Pro
ip velx 1, eptxx 0, while the old GNU gave ip velx 0.5, eptxx -0.01.

## Implementation

- `top.cc` `top()`, immediately after the equilibrium `for (iteration...)`
  loop and before the old hoisted stress-only block:
  `iteration_start(); material_state_post_pass = 1; element_loop();
  material_state_post_pass = 0; parallel_sys_routine(
  &parallel_new_dof_diagonal ); mpc_node_apply();`
  It runs when `nuknwn>0`, `options_solver != -DIAGONAL`,
  `options_elementloop == -YES`, `!find_local_softvar`, no
  `options_nonlocal_softvar`, and `post_solve_state_pass_gated(icontrol)`
  returns 0. Otherwise the hoisted stress-only sync
  (`nodal_stress_post_solve_sync`) runs exactly as before.
- `top.cc` global `long int material_state_post_pass` (declared in
  `tochnog.h` and `tochnog-mod.h`). `elem.cc` `elem()` skips the sparse
  ELEMENT_MATRIX_VALUES / GROUP_MATRIX_VALUES store while it is set: no
  solve follows the pass. The stored matrices stay those of the last
  solve, which is also what a step-close `control_eigen` uses.
- What the pass keeps: the full element loop (material state, element
  force records, NODE_RHSIDE/NODE_LHSIDE, i.e. the residual at the final
  iterate used for the reported reactions), the non-principal nodal update
  (`parallel_new_dof_diagonal`: absolute stress recovery, incremental
  correction of the other non-principal dofs; the correction form
  `rate - (new_dof - past)/dt` makes the second application idempotent)
  and the mpc re-sync. `post_node_rhside_fixed_free()` is not re-run:
  POST_NODE_RHSIDE_* keep the values of the last iteration.
- The element tangent is still computed in the element scratch arrays (it
  is a by-product of the constitutive routines); only its store is skipped.

### Gate: `post_solve_state_pass_gated()` (top.cc)

Returns 1 (pass OFF) for:

- slide: NODE_SLIDE or any record of data class SLIDE;
- contact: any record of data class CONTACT, unless CONTACT_APPLY or
  CONTROL_CONTACT_APPLY (this icontrol) is -no (same test as
  `contact.cc` `parallel_contact()`).

Reason (prototype A/B without the gate): Pro `dynamic3` (truss, prescribed
acceleration) and `slide4` (slide friction) moved AWAY from the
Professional. Contact and slide forces are added outside the element loop
(`parallel_contact()`, `slide()`), and the inertia terms use the
velocity/acceleration update of the iteration, so the pass does not yet
update that physics consistently.

**Widening criterion** (project owner 2026-10-10): a gated family is
enabled only after the pass is extended to its physics (for example
calling `slide()` / `parallel_contact` consistently, or a consistent
inertia update) AND the probes plus the full corpus A/B on omen show no
test moving away from the Professional. Investigations: BACKLOG C52
(`dynamic3`, inertia/truss) and C53 (`slide4`, friction state).

**Inertia widened (F4 T6d, BACKLOG C52, 2026-10-11)**: the inertia test
(`materi_acceleration`, OPTIONS_INERTIA, CONTROL_OPTIONS_INERTIA) was
removed from the gate. The two blockers were not the pass: Pro `dynamic3`
was the truss memory default / sign defect C11 (= Professional with and
without the gate), and the test-2014 `iface_gen`, `iface_gen_geom`,
`iface_gen_method` break of the first attempt (node 6 vely -1.70554, rc 1)
was the interface Mohr-Coulomb elastic-trial history (fixed in
`interface.cc`, `manual-developer/group_interface.md` record 2026-10-11).
With both fixes these tests are rc 0 with the pass. An interface +
`inertia_apply` probe keeps an older GNU/Professional gap that
the pass does not cause (BACKLOG C61). Slide and contact stay gated (C53).

## Verification

- Probes (Professional 02-08-2026, scratch `f4t5/probes`): 1D stretch
  velocity-only 1 and 6 steps, with `materi_displacement`, two plastic
  bars 1 and 5 steps, `spring4` rigid translation. All are equal to the
  Professional after the change, e.g. plastic bars node 1 eppxx 0.01 /
  sigxx 0.015 (Pro 0.0100001 / 0.0149999; old 0 / 0).
- Matrix store vs no matrix store: the TN_C35 prototype (full store) and
  the final build give byte-identical `.dbs` on the probes and on every
  non-gated corpus test (742), except 2-thread round-off of
  `bicg_csr_threads1`.
- Tests: `test-2026/post_solve_state_stretch1` (node 2 velx 0.5 in all six
  steps, `velx2.his`), `post_solve_state_plastic1`,
  `post_solve_state_rigid1` (RED on the old code, rc 1).
- Corpus A/B: see SEGUIMIENTO (F4 T5 row) and `odd/tasks/f4-investigations.md`.

## Hardcoded parameters / pending refactorings

- The gate is a hard-coded list of data classes and records. When a family
  is widened, remove its test from `post_solve_state_pass_gated()`.
- The pass and the stress-only sync coexist. Once every family is covered,
  the stress-only block (and `stress_recovery_post_pass`) can be removed.
- Cost: one extra element loop per step (no matrix store, no solve).
  Models with many cheap steps and few iterations pay the most (Pro
  `ground11_nonsaturated` 84 -> 105 s alone).

## Differences with the Professional version

- Gated families (slide, contact since F4 T6d) keep the one-iteration lag
  (DIVERGENCES class B `POSTPASS-GATE`).
