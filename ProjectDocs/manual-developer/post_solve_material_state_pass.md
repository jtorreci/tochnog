# post_solve_material_state_pass

Status: **IMPLEMENTED 2026-10-10** (F4 T5, BACKLOG C35, audit K02;
project-owner decision 2026-10-10). This page replaces
`nodal_stress_post_solve_sync` for every model (since F4 T6e, 2026-10-11,
no physics is gated out).

## Problem

`top.cc` `top()` iterates `element_loop(); solve(); parallel_new_dof_diagonal`.
The element loop of iteration k evaluates the material on the PRE-solve
iterate (velocity of iteration k-1). It writes the integration-point state
(ELEMENT_INTPNT_DOF VERSION_NEW: stress, strains, plastic strain, kappa, history
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
  material_state_post_pass = 0; parallel_sys_routine( &parallel_contact );
  slide(); parallel_sys_routine( &parallel_new_dof_diagonal );
  mpc_node_apply();` (contact and slide re-applied since F4 T6e C53)
  It runs when `nuknwn>0`, `options_solver != -DIAGONAL`,
  `options_elementloop == -YES`, `!find_local_softvar`, no
  `options_nonlocal_softvar`. Otherwise the hoisted stress-only sync
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
  and the mpc re-sync. `post_node_rhside_fixed_free()` is re-run after
  the pass (C110, 2026-11-05): POST_NODE_RHSIDE_* are those of the final
  residual.
- The element tangent is still computed in the element scratch arrays (it
  is a by-product of the constitutive routines); only its store is skipped.

### Former gate: `post_solve_state_pass_gated()` (top.cc, removed F4 T6e)

Until 2026-10-11 it returned 1 (pass OFF) for:

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

**Slide and contact widened, gate removed (F4 T6e, BACKLOG C53,
2026-10-11)**: the pass now calls `parallel_sys_routine( &parallel_contact )`
and `slide()` after its element loop, in the order of an equilibrium
iteration (the iteration also runs `parallel_new_dof_before` between them;
the pass does not, see below). `post_solve_state_pass_gated()` and its
declarations in `tochnog.h` / `tochnog-mod.h` were deleted. The Pro
`slide4` gap of the first attempt needed, besides the slide/contact
re-application, two slide-law fixes (`slide.cc`, commit `fix(slide)`
of the same task): an absent `slide_plasti_tension` defaulted to
`sig_t = 0` (a compression-only gap) and the friction cap was
`max(c + mu*S, 0)` instead of the Professional `c + mu*|S|`. With all of
it Pro `slide4` = Professional (node 1 `node_slide_force`
-2.11291867e-3 / -1.21989198e-3 vs Pro -2.11291887e-3 / -1.21989209e-3;
HEAD rc 1, -2.0602e-3), the cohesion variant `slide_friction_cohesion1`
-2.06548203e-3 (Pro -2.06548222e-3) and the Pro-syntax contact_block
`contact_block_pass1` node_dof within 7e-13 of the Professional. The
read-only investigation A/B-tested the patch on omen (`f4t6e_runs`). With
the final binary the slide/contact focus set (Pro slide1-4, contac1/2/4,
test-2014 contact*, node_slide, slide_axi, test-2026 contact_*) is 18/18
rc 0 (HEAD: slide4 rc 1); build_safe 216 OK / 0 FALLO, runs 226/226.

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

- C62 (F4 T6f LOW queue Q3, 2026-10-17, resolved): the pass now runs
  the same nodal sequence as an equilibrium iteration. `top.cc` keeps a
  file-static snapshot of `NODE_RHSIDE` taken right after `bounda()` in
  every iteration (`bounda_rhside_snapshot_take()`: the `bounda_force`
  loads, which `bounda()` ASSIGNS before the element loop) and restores
  it after the pass's `iteration_start()` (`_restore()`, skipped when the
  node count changed); the pass then calls `parallel_new_dof_before`
  between `parallel_contact` and `slide()` (node_force, force_point,
  node_mass/damping/stiffness, nodal gravity, `NODE_INERTIA` on the
  final iterate). Nothing is applied twice: both only touch the
  velocity slots of `NODE_RHSIDE`/`NODE_LHSIDE`, and for a non-diagonal
  solver `parallel_new_dof_diagonal` updates only the non-principal dofs
  (the `-diagonal` solver never takes the pass). `bounda()` itself is not
  re-run (it re-prescribes dofs, `bounda_until_force` scales velocities).
  Tests (Professional targets, RED on `feafe50`):
  `test-2026/node_rhside_post_pass_node_force1`,
  `_node_mass1`, `_bounda_force1`.
- The pass and the stress-only sync coexist. Once every family is covered,
  the stress-only block (and `stress_recovery_post_pass`) can be removed.
- Cost: one extra element loop per step (no matrix store, no solve).
  Models with many cheap steps and few iterations pay the most (Pro
  `ground11_nonsaturated` 84 -> 105 s alone).

## Differences with the Professional version

- None from the pass since F4 T6e (DIVERGENCES `POSTPASS-GATE` removed).
  C62 resolved 2026-10-17. C63 resolved 2026-11-05 by C110
  (`post_node_rhside_fixed_free()` re-run after the pass; Pro `slide4`
  = Professional).
