# control_restart

## Where

- `restart.cc`, `restart()`: called from `top.cc` (control-index start,
  after `refine_locally()` and before `unknown_reset()`).
- `restart.cc`, `restart_reset_element_state()` (static, F4 T6c).

## Implementation

For `control_restart -yes` at the current control index, `restart()`:

1. copies `time_current` from `VERSION_START` to `VERSION_NORMAL`;
2. `db_copy( NODE_START_REFINED, NODE )` and
   `db_copy( NODE_DOF_START_REFINED, NODE_DOF )` (start coordinates and
   dofs of the current refined mesh, written in `top.cc`);
3. `restart_reset_element_state()`: sets `ELEMENT_DOF_INITIALISED` to 0
   for every active element in `VERSION_NORMAL` and, when the version
   exists, `VERSION_NEW`.

With `options_element_dof -yes` (default, `initia.cc`) `elem.cc`
restores the remembered integration-point stress/strain/history block
from `ELEMENT_DOF` only when `element_dof_initialised` is set; with the
flag cleared the IP values are interpolated from the (reset) nodal dofs,
which gives a fresh integration-point state. `top.cc` sets the flag back
to 1 after the first accepted step. This mirrors the Professional, which
deletes `element_dof` / `element_intpnt_dof` on a restart.

Verification (F4 T6c): `validation-suite/test-2026/control_restart_reset_ip1`
(Drucker-Prager) and `_ip2` (elastic), targets = Professional 02-08-2026;
RED on the old binary (rc 1, sigyy sign flipped), GREEN after the fix.

## Element state not covered

The flag only resets the integration-point dofs. Other element-level
history records are untouched:

- `ELEMENT_TRUSS_FORCE`: kept, = Professional (probe 1.5 in both).
- `ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE`: result = Professional
  (probe -9.75297 in both).
- Interface (`ELEMENT_INTERFACE_STRAIN_NORMAL`, `_FORCE_NORM`,
  `_FORCE_TANG`, `_INTPNT_STRESS`, tension status) and spring
  (`ELEMENT_SPRING_FORCE`, `-total_linear` and `-updated`): kept, =
  Professional (F4 T6f LOW queue Q6, 2026-10-17, probes
  `scratchpad/lowq/q6/`: the Professional keeps the accumulated history
  too: interface -1.5 / 1.5 after a restart, fresh run -0.5 / 0.5;
  spring 0.15 vs 0.05). Tests `test-2026/control_restart_interface_history1`,
  `control_restart_spring_history1` (pass before and after: they pin
  the behaviour). Side findings (not restart): the Professional writes
  an `element_spring_strain` record the GNU lacks, and the GNU
  `element_interface_intpnt_strain` holds the last step's increment
  (SEGUIMIENTO PENDIENTE).
- Not verified (deferred, class B C58): `ELEMENT_CONTACTSPRING_FORCE`,
  `ELEMENT_TENDON_STRAIN/_STRESS`,
  `ELEMENT_INTPNT_MATERI_STRESS_EFFECTIVE` (damage).
- `number_timesteps`: the Professional resets it; the GNU has no such
  record.

## Pending refactorings

- `ELEMENT_DOF_INITIALISED` is a 0/1 integer flag per element; the reset
  relies on `elem.cc` interpolation semantics rather than an explicit
  "fresh state" API.
