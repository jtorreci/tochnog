# post_apply — developer

## Implementation

- **Enums**: `POST_APPLY` (existing), `CONTROL_POST`, `CONTROL_POST_APPLY`,
  `CONTROL_POST_ELEMENT_FORCE` (`tochnog.h` / `tochnog-mod.h`, after
  `POST_APPLY`); records INTEGER length 1, data_class CONTROL
  (`database.cc`).
- **Gate**: `step_close()` (`top.cc`): `post_off = materi_apply_gate_off(
  POST_APPLY, CONTROL_POST_APPLY )` (control at `ICONTROL` overrules the
  global record, default -yes) skips `post( task )` and `calculate()`
  (`post_calcul`). The timestep path skips `post_element_force_calculate()`
  under the same gate and under `control_materi_gate_off(
  CONTROL_POST_ELEMENT_FORCE )`. `post_node_rhside_fixed_free()` (in the
  iteration loop) is not gated.
- **Stale flags** (globals `post_stale`, `post_element_force_stale`,
  `tochnog.h` / `tochnog-mod.h`, defined in `initia.cc`): set when the last
  evaluation was skipped. A control index without timesteps recomputes a
  stale `post_element_force` (the `post()` call of its `step_close( NO )`
  already runs); `exit_tn( -YES )` (`miscel.cc`) redoes stale records
  (`post( NO )`, `calculate()`, `post_element_force_calculate()`) before
  the final database print, unless the global `post_apply` is `-no`. The
  Professional re-evaluates at the end of every run; redoing only the
  stale records gives the same output and leaves runs without these
  switches byte-identical.

## Verification (Professional binary 02-08-2026 as oracle)

Quad4 (E 1000, nu 0.3), top vely ramp 0 -> 0.02 over t 0..1, two control
indices of 5 steps, `post_point 0 0.5 0.5`, `post_element_force` on the
bottom edge; series from `control_print_history` per index. Equal series in
both binaries for: no switch; `post_apply -no` (no record, no history);
`control_post_apply 0 -no` (history starts at t 0.6); `control_post_apply 1
-no` (disy frozen 0.0015, force 3.2967); global -no + control 1 -yes;
`control_post_element_force 0/1 -no` (only the force frozen);
`control_post 2 -yes` / `-no` / absent with an index 2 without timesteps
(index 2 evaluates: 0.0055, 12.0879); `control_post 1 -no` (no effect).
Final `.dbs`: 0.0055 / 12.0879 in every case except `post_apply -no`.
Differential vs `4b80565` on the 212 test inputs with post records: only
the new tests and `post_apply1` change.

## Pending

- `post_integrate` runs only with task YES; `control_post -yes` does not
  force it on an index without timesteps (no time increment there).
