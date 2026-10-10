# post_solve_material_state_pass

Status: **IMPLEMENTED 2026-10-10** (F4 T5, BACKLOG C35, audit K02). No
keyword: it is part of the iteration scheme and always on, except for the
model families listed under "Gate".

## What it does

Within a time step Tochnog repeats equilibrium iterations: evaluate the
material on the current velocity, assemble, solve for a new velocity. The
material evaluation of an iteration uses the velocity from BEFORE that
iteration's solve. Until 2026-10-10 the state written at the end of the
step was therefore one iteration behind the solved velocity. This applied to
the integration-point stress, strain, plastic strain and hardening, the
nodal strain/plastic strain/kappa dofs and the next step's history. With
few iterations the reported state could be visibly wrong. For example, a
rigid translation gave strain -0.01, and a one-iteration plastic bar did not
yield.

After the last equilibrium iteration of every step, the material is now
evaluated once more on the final (solved) velocity. That state is the
output of the step and the history of the next one. This is what the
Professional does. The extra pass does not solve and does not build a
system matrix, so it costs about one extra element loop per step.

## Why it is useful

- Results with a small number of equilibrium iterations
  (`control_timestep_iterations 1`, explicit-like runs) now report the
  state of the solved velocity, as in the Professional.
- The reported velocity of free nodes and the strain/stress of the same step
  are now consistent.

## Input

None. The pass is automatic.

## Gate (temporary difference with the Professional)

For these families the pass is OFF and the previous scheme is kept
unchanged (only the nodal stress is re-synchronised, see
`nodal_stress_post_solve_sync`):

- inertia: `materi_acceleration`, `inertia_apply -yes` /
  `options_inertia -yes`, or any `control_inertia_apply -yes`;
- slide: any `slide_*` record or `node_slide`;
- contact: any `contact_*` record, unless `contact_apply -no` /
  `control_contact_apply -no` for the step;
- also the `-diagonal` solver, `options_elementloop -no` and the nonlocal
  softening variable.

## Example

```
number_of_space_dimensions 1
materi_velocity
materi_velocity_integrated
materi_strain_total
materi_stress
end_initia
node 1 0
node 2 1
node 3 2
element 1 -bar2 1 2
element 2 -bar2 2 3
bounda_dof 0 1 -velx
bounda_time 0 0.
bounda_dof 1 3 -velx
bounda_time 1 1.
group_type 0 -materi
group_materi_elasti_young 0 1.
control_timestep 0 0.01 0.06
control_timestep_iterations 0 1
end_data
```

Node 2 now has `velx` 0.5 in every step and `eptxx = sigxx = 0.03` at the
end, the same as the Professional 02-08-2026. Before, node 2 velx went
0.5, 1, 1, 0.5, 0, 0. See `test-2026/post_solve_state_stretch1`,
`post_solve_state_plastic1` and `post_solve_state_rigid1`.

## Differences with the Professional version

- Gated families (inertia, slide, contact): the GNU keeps the previous
  scheme, so with few iterations their state still lags one iteration.
  This is a temporary divergence (DIVERGENCES class B, `POSTPASS-GATE`). The
  measured cases are Pro `dynamic3` (truss with prescribed acceleration) and
  `slide4` (slide friction state): there, the ungated pass moved the result
  away from the Professional. A family is enabled only after the pass
  covers its physics and the full corpus shows no test moving away from the
  Professional.
