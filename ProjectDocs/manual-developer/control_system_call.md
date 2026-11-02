# control_system_call — developer

## Implementation

- Record `CONTROL_SYSTEM_CALL` (INTEGER, data_length 1, data_class CONTROL;
  `database.cc` after `CONTROL_SUPPORT_EDGE_NORMAL_STIFFNESS_FREEZE`, enum in
  `tochnog.h` / `tochnog-mod.h`).
- `control_system_call( icontrol )` in `top.cc` (end of file). Called from
  `top()` right before `step_close( YES, ... )` of a converged timestep and
  right before `inverse_calculation( ..., INVERSE_DETERMINE_SENSITIVITY )`
  at the end of the control index (both the timestep and the no-timestep
  branch).
- The script is started with `posix_spawn()` (no shell, no `system()`): the
  input supplies only the integer, formatted with `%ld`. File action: fd 2
  opened on `/dev/null`. On `ENOEXEC` (no `#!` line) it retries with
  `/bin/sh ./tochnog_system_call.sh N`. `waitpid()` waits for the child
  (EINTR retried); the exit status is ignored. `cout` and `stdout` are
  flushed before the spawn so the script output keeps its order.
- `access( X_OK )` check first; a static flag limits the missing-script
  warning to one per run.
- Negative values are already rejected by the input parser (`-3` is read
  as an unknown word: "I do not know -3"), as the Professional.

## Hardcoded parameters

- Script name `./tochnog_system_call.sh` (manual), stderr target `/dev/null`.

## Verification

Professional 02-08-2026 probes (quad4 `control_post1` geometry, indices 0
(3 steps), 1 (2 steps), 2 (no timesteps), script logging `$1` and the
history lines already written): 7 x4 (his 0,1,2,3), 8 x3 (3,4,5), 9 x1 (5);
stdout kept, stderr dropped, exit 3 ignored; missing / non-executable
script silent; `-3` rejected. The Professional only finds the script with
`.` in the `PATH` (`system("tochnog_system_call.sh N")`). GNU log identical.

## Pending

- `control_repeat`: the end-of-index call runs once after the repeated
  index has finished (the repeat jumps back before it); not measured on the
  Professional.
- MS Windows `.bat` variant not supported.
